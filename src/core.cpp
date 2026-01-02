//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#include "core.hpp"


core::core()
{

    _start_ts = nullptr;
    _end_ts = nullptr;
    _interpolation_method = interp_alg::tpspline;
    radius = 0;
    N=0;

    point_mode.enable = false;
    point_mode.use_specific_station = false;
    point_mode.forcing = "";

    //default logging level
    _log_level = debug;

    //chkpoointing
    _use_netcdf=false;

    _metdata= nullptr;

    clean_exit = true;

}

core::~core()
{
    // clean up all the modules to ensure that Tpetra:~Map() is called prior to MPI_Finalize as per GitHib Issue #2372
    for(auto& itr : _chunked_modules)
    {
        for(auto& jtr : itr )
        {
            jtr.reset();
        }
    }
    for(auto& itr : _modules)
    {
        itr.first.reset();
    }


    // flush the logs and close the fhandles
    spdlog::shutdown();

    try
    {
        auto job_name = _hpc_scheduler_info.job_name;

        const auto path = output_folder_path / ("logs_" + job_name);
        try{ boost::filesystem::create_directories(path); } catch (...) { /* already exists */ }

        boost::filesystem::rename(log_file_path, path / log_file_path.filename());
    }
    catch(const std::exception& e)
    {
        // if this copying goes wrong we are super out of options
        std::cout << e.what() << std::endl;
    }

}

void core::config_options( pt::ptree &value)
{
    SPDLOG_DEBUG("Found options section");

    _find_and_insert_subjson(value);

    //setup debugging level, default to debug if not specified
    std::string s = value.get("debug_level", "debug");

    if (s == "debug")
        _log_level = debug;
    else if (s == "warning")
        _log_level = warning;
    else if (s == "error")
        _log_level = error;
    else if (s == "verbose")
        _log_level = verbose;

    SPDLOG_DEBUG("Setting log severity to {}", s);


    // guard against deprecated config files
    boost::optional<std::string> interp = value.get_optional<std::string>("interpolant");
    boost::optional<std::string> ssr = value.get_optional<std::string>("station_search_radius");
    boost::optional<std::string> sNn = value.get_optional<std::string>("station_N_nearest");
    if(interp || ssr || sNn)
    {
        SPDLOG_ERROR("Interpolation and N station config options have moved to the forcing section. Please update "
                     "your configuration file.");
        CHM_THROW_EXCEPTION(config_error, "Old style config");
    }


    // custom start time
    boost::optional<std::string> start = value.get_optional<std::string>("startdate");
    if (start)
    {
        _start_ts = new boost::posix_time::ptime(boost::posix_time::from_iso_string(*start));
        SPDLOG_DEBUG("User-specified startdate: {}", boost::posix_time::to_simple_string(*_start_ts));
    }

    // custom start time
    boost::optional<std::string> end = value.get_optional<std::string>("enddate");
    if (end)
    {
        _end_ts = new boost::posix_time::ptime(boost::posix_time::from_iso_string(*end));
        SPDLOG_DEBUG("User-specified endate: {}",  boost::posix_time::to_simple_string(*_end_ts));
    }


    // point mode options
    try
    {
        auto pm = value.get_child("point_mode");

        point_mode.enable = true;
        point_mode.use_specific_station = false;

        try
        {
            // if we ask for a specific forcing point, then /only/use that point.
            // Only works with ascii mode so we will need to check for this later
            point_mode.forcing = pm.get<std::string>("forcing");
            point_mode.use_specific_station = true;
            if (_interpolation_method != interp_alg::nearest_sta)
            {
                _interpolation_method = interp_alg::nearest_sta;
                SPDLOG_WARN( "Station select has been changed to nearest station because a single point mode station was requested");
            }
        }
        catch (pt::ptree_bad_path& e)
        {
            //pass
        }

        _global->_is_point_mode = true;

    }
    catch(pt::ptree_bad_path &e)
    {
        point_mode.enable = false; // we don't have point_mode
    }

    auto notify_sh = value.get_optional<std::string>("notification_script");
    if(notify_sh)
    {
        _notification_script = *notify_sh;
    }

}

void core::config_module_overrides( pt::ptree &value)
{
    _find_and_insert_subjson(value);
    SPDLOG_DEBUG("Found dependency override section");
    for (auto &itr : value)
    {
        std::string A = itr.first.data();
        std::string B = itr.second.data();

        SPDLOG_WARN("Removing depency of {} from {}", B, A);
        _overrides.push_back(std::make_pair(A, B));

    }
}

void core::config_modules(pt::ptree &value, const pt::ptree &config, std::vector<std::string> remove,
                          std::vector<std::string> add)
{
    SPDLOG_DEBUG("Found modules section");
    int modnum = 0;
    //loop over the list of requested modules
    // these are in the format "type":"ID"
    std::set<std::string> modules; //use a set to avoid duplicates (ie, user error!)


    for (auto &itr : value)
    {
        std::string module = itr.second.data();

        //see if the module is one in the remove list
        //if it isn't, add it
        if (std::find(remove.begin(), remove.end(), module) == remove.end())
        {
            modules.insert(module);
        }
        else
        {
            SPDLOG_DEBUG("Removed module {}", module);
        }

    }

    //straight up add anything
    for (auto &itr : add)
    {
        modules.insert(itr);
        SPDLOG_DEBUG("Inserted module {} from cmdl", itr);
    }


    for (auto &itr : modules)
    {

        std::string module_name = itr;
        SPDLOG_DEBUG("Module ID={}", module_name);


        //try grabbing a config for this module, empty string default
        pt::ptree cfg;
        try
        {
            cfg = config.get_child(module_name);
        } catch (pt::ptree_bad_path &e)
        {
            SPDLOG_DEBUG("No config for {}", module_name);
        }

        boost::shared_ptr<module_base> module = module_factory::create(module_name, cfg);
        //internal tracking of module initialization order
        module->IDnum = modnum;

        //assign the internal global param pointer to our global
        module->global_param = _global;

        //get the parameters that this module will provide
        // we need this now before the mesh call as in the mesh call we will build the static mph hashtable for all parameters
        for(auto& p: *(module->provides_parameter())) {
            _provided_parameters.insert(p);
        }

        for(auto& p: *(module->provides_vector())) {
            _provided_var_vector.insert(p);
        }


        modnum++;
        _modules.push_back(
                std::make_pair(module, 1)); //default to 1 for make ordering, we will set it later in determine_module_dep
    }


    if (modnum == 0)
    {
        CHM_THROW_EXCEPTION(no_modules_defined, "No modules defined. Aborting");
    }
}

void core::config_checkpoint( pt::ptree& value)
{
    SPDLOG_DEBUG("Found checkpoint section");

    _checkpoint_opts.do_checkpoint = value.get("save_checkpoint",false);

    //this uses the output path defined in config_output, so this must run after that.
    if(_checkpoint_opts.do_checkpoint)
    {

        auto dir = "checkpoint";

        // Create empty folder
        _checkpoint_opts.ckpt_path = output_folder_path / dir;
        boost::filesystem::create_directories(_checkpoint_opts.ckpt_path);

        // check for the old naming and bail
        auto tmp = value.get_optional<size_t>("frequency");
        if(tmp)
        {
            CHM_THROW_EXCEPTION(chm_error, "checkpointing.frequency has been renamed to checkpointing.on_frequency");
        }

        _checkpoint_opts.frequency = value.get_optional<size_t>("on_frequency");
        if (_checkpoint_opts.frequency)
        {
            SPDLOG_DEBUG("Checkpointing every {} timesteps" , *(_checkpoint_opts.frequency));
        }

        _checkpoint_opts.on_last = value.get_optional<bool>("on_last");
        if (_checkpoint_opts.on_last && *(_checkpoint_opts.on_last))
        {
            SPDLOG_DEBUG("Checkpointing on last timestep");
        }

        // checkpoint before we reach a wall clock limit
        _checkpoint_opts.on_outta_time = value.get_optional<bool>("on_wallclock_limit");
        if(_checkpoint_opts.on_outta_time && *(_checkpoint_opts.on_outta_time))
        {
            auto tmp = value.get_optional<size_t>("minutes_of_wallclock");

            if(tmp)
            {
                _checkpoint_opts.abort_when_wallclock_left = boost::posix_time::minutes(*tmp);
            }

            SPDLOG_DEBUG("Checkpointing if < {} of wallclock left",
                         boost::posix_time::to_simple_string(_checkpoint_opts.abort_when_wallclock_left));

            if(!_hpc_scheduler_info.has_wallclock_limit)
            {
                SPDLOG_ERROR("Wallclock limit checkpointing is enabled, but no CHM_WALLCLOCK_LIMIT "
                             "environment variable has been detected.");
                CHM_THROW_EXCEPTION(chm_error, "Missing wallclock envar required for checkpoint");
            }
        }

        auto specific_datetime = value.get_optional<std::string>("specific_datetime");
        if(specific_datetime)
        {
            _checkpoint_opts.specific_datetime = boost::posix_time::ptime(boost::posix_time::from_iso_string(*specific_datetime));
            SPDLOG_DEBUG("Checkpointing every {} timesteps" , boost::posix_time::to_simple_string(*(_checkpoint_opts.specific_datetime)));
        }

        auto specific_time = value.get_optional<std::string>("specific_time");
        if(specific_time)
        {
            // ptime needs a real date, but we will only ever be checking the minute and hour
            auto time = "3000-01-01 " + *specific_time + ":00";
            _checkpoint_opts.specific_time = boost::posix_time::ptime(boost::posix_time::time_from_string(time));
            SPDLOG_DEBUG("Checkpointing every {} timesteps" , boost::posix_time::to_simple_string(*(_checkpoint_opts.specific_time)));
        }

        if( !_checkpoint_opts.on_last &&
            !_checkpoint_opts.frequency &&
            !_checkpoint_opts.on_outta_time &&
            !_checkpoint_opts.specific_datetime &&
            !_checkpoint_opts.specific_time)
        {
            SPDLOG_ERROR("Checkpointing is enabled but no `on_*` selection has been specified.\n"
                         "Please see https://chm.readthedocs.io/en/develop/configuration.html#checkpoint\n");
            CHM_THROW_EXCEPTION(config_error, "Missing checkpoint options");
        }


    }

    // automatically resume from the most-recent (simulation time) checkpoint file that exists in the output directory
    auto auto_resume = value.get_optional<bool>("auto_resume");

    // file to resume off of
    auto file = value.get_optional<std::string>("load_checkpoint_path");

    //before spending any effort on this make sure we don't have contradictory options
    if(auto_resume && *auto_resume && file)
    {
        SPDLOG_ERROR("Checkpoint options auto_resume and file cannot both be specified.");
        CHM_THROW_EXCEPTION(config_error, "Error in checkpoint config");
    }

    if(*auto_resume)
    {

        // match the filename of checkpoint_20170901T070000.np1.json
        // with the iso date as the single capture
        boost::regex file_pattern("checkpoint_(\\d{8}T\\d{6})\\.np\\d+\\.json");

        boost::posix_time::ptime newest_time;
        boost::filesystem::path newest_file; // will end up with the most recent chkpt to resume from
        boost::filesystem::path const dir(output_folder_path  / "checkpoint" );

        for(const auto& entry : boost::filesystem::directory_iterator(dir))
        {
            if (boost::filesystem::is_regular_file(entry.status()))
            {
                auto filename = entry.path().filename().string();
                boost::smatch match;

                if (boost::regex_match(filename, match, file_pattern))
                {
                    std::string datetime_str = match[1];
                    auto file_time = boost::posix_time::from_iso_string(datetime_str);

                    // keep looking until we find the most recent date
                    if (newest_file.empty() || file_time > newest_time)
                    {
                        newest_time = file_time;
                        newest_file = entry.path();
                    }
                }
            }
        }

        if(newest_file.empty())
        {
            SPDLOG_WARN("Could not find a checkpoint to resume from.");
        }
        else
        {
            file = newest_file.string();
        }

    }



    if (file)
    {
        _checkpoint_opts.load_from_checkpoint = true;

        boost::filesystem::path ckpt_path = *file;
//        ckpt_path = boost::filesystem::canonical(ckpt_path);

        auto chkp = read_json(ckpt_path.string());

        size_t csz = 1;
        size_t rank = 0;

        csz = _comm_world.size();
        rank = _comm_world.rank();


        if( csz != chkp.get<size_t>("ranks") )
        {
            CHM_THROW_EXCEPTION(config_error, "Checkpoint file was saved with a different number of ranks");
        }

        try
        {
            _global->timestep_counter = chkp.get<size_t>("timestep_counter");
        }catch (...)
        {
            CHM_THROW_EXCEPTION(config_error,
                "This checkpoint was saved with an older version and does not have the field timestep_counter");
        }

        boost::filesystem::path ckpt_nc_path;
        try
        {
          size_t i = 0;

            for(auto &itr : chkp.get_child("files"))
            {
                if( i == rank)
                {
                    ckpt_nc_path = itr.second.data();
                    break;
                }

                ++i;
            }
        }
        catch(pt::ptree_bad_path &e)
        {
          CHM_THROW_EXCEPTION(config_error, "Error reading list of checkpoint files");
        }

        ckpt_nc_path =  ckpt_path.parent_path() / ckpt_nc_path;
        SPDLOG_DEBUG("Rank {} using checkpoint restore file {}", rank, ckpt_nc_path.string());
        _checkpoint_opts.in_savestate.open(ckpt_nc_path.string());
    }
}
void core::config_forcing(pt::ptree &value)
{
    SPDLOG_DEBUG("Found forcing section");

    std::string ia = value.get<std::string>("interpolant","spline");

    if( ia == "spline")
    {
        _interpolation_method = interp_alg::tpspline;
    }
    else if (ia == "idw")
    {
        _interpolation_method = interp_alg::idw;
    }
    else if (ia == "nearest")
    {
        _interpolation_method = interp_alg::nearest_sta;
    }
    else
    {
        SPDLOG_WARN("Unknown interpolant selected, defaulting to spline");
    }

    if(value.get_optional<double>("station_N_nearest"))
    {
        CHM_THROW_EXCEPTION(config_error, "station_N_nearest option is renamed num_stations_to_use and moved to the forcing section.");
    }

    auto radius = value.get_optional<double>("station_search_radius");
    auto N = value.get_optional<int>("num_stations_to_use");

    if(radius && N)
    {
        CHM_THROW_EXCEPTION(config_error, "Cannot have both station_search_radius and num_stations_to_use set.");
    }

    if(!radius)
    {
        if(!N)
        {
            N = 5;
        }

        if (_interpolation_method == interp_alg::nearest_sta)
        {
            N = 1;
            SPDLOG_DEBUG("Setting N=1 for nearest interpolant");
        }

        if( (*N < 2) && (_interpolation_method != interp_alg::nearest_sta)) // Required more than 1 station if using spline or idw
        {
            CHM_THROW_EXCEPTION(config_error, "num_stations_to_use must be >= 2 if spline or idw is used. N = " + std::to_string(*N));
        }

        SPDLOG_DEBUG("Using N={} nearest stations as default.", *N);
    }

    // This needs to be initialized with the mesh prior to the forcing being loaded because
    // metdata needs to know about the meshes' coordinate system and bounding box to only load nearby stations
    _metdata = std::make_shared<metdata>(_mesh, output_folder_path, N);

    // Set get_stations to be the correct lookup function
    if(radius)
    {
        _metdata->get_stations = boost::bind( &metdata::get_stations_in_radius,_metdata,boost::placeholders::_1,boost::placeholders::_2, *radius);
    }
    else
    {
        _metdata->get_stations = boost::bind( &metdata::nearest_station,_metdata,boost::placeholders::_1,boost::placeholders::_2, *N);
    }

    SPDLOG_DEBUG("Reading meta data from config");



    //positive offset going west. So the normal UTC-6 would be UTC_offset:6
    _global->_utc_offset = value.get("UTC_offset",0);
    SPDLOG_DEBUG("Applying UTC offset to ALL forcing files. UTC+{}", std::to_string(_global->_utc_offset));

    //need to determine if we have been given a netcdf file
    _use_netcdf = value.get("use_netcdf",false);


    timer c;
    c.tic();
    size_t nstations = 0;
    //we need to treat this very differently than the txt files

    if(_use_netcdf)
    {
        std::string file = value.get<std::string>("file");
        std::map<std::string, boost::shared_ptr<filter_base> > netcdf_filters;
        try
        {
            auto filter_section = value.get_child("filter");

            for (auto &jtr: filter_section)
            {
                auto filter_name = jtr.first.data();
                auto cfg  = jtr.second;

                boost::shared_ptr<filter_base> filter = filter_factory::create(filter_name,cfg);
                filter->init();
                netcdf_filters[filter_name] = filter;
            }

        } catch (pt::ptree_bad_path &e)
        {
            //ignore bad path, means we don't have a filter
        }

        // Check if we are loading a list from a json file of netcdffiles
        if(file.find(".json") != std::string::npos)
        {
            _metdata->load_from_listof_netcdf(file,netcdf_filters);
        }
        else
        {
            // this delegates all filter responsibility to metdata from now on
            _metdata->load_from_netcdf(file, netcdf_filters);
        }

        nstations = _metdata->nstations();

        // we need to estimate elevations of the forcing assuming it's surface level from the mesh
        if(_metdata->missing_z())
        {
            SPDLOG_WARN("No geopotential height field found in the netcdf file. Using the height of nearest triangle to estimate. This is almost certainly NOT what you want");

            for(size_t i = 0; i < _mesh->size_local_faces(); i++)
            {
                for (auto face = _mesh->face(i); auto& s : face->stations())
                {
                    if(s->z() == -9999)
                    {
                        s->z(face->get_z());
                    }
                }

            }
        }
    } else
    {
        std::vector<metdata::ascii_metdata> ascii_data;

        for (auto &itr : value)
        {
            if(itr.first != "UTC_offset" &&
                itr.first != "num_stations_to_use" &&
                itr.first != "interpolant")
            {
                metdata::ascii_metdata data;

                std::string station_name = itr.first.data();
                auto& station = value.get_child(station_name);

                data.id = station_name;

                try
                {
                    data.longitude = station.get<double>("longitude");
                    data.latitude = station.get<double>("latitude");
                }catch(boost::property_tree::ptree_bad_path& e)
                {
                    CHM_THROW_EXCEPTION(forcing_error,"Station " + station_name + " missing lat/long coordinate.");
                }

                data.elevation = station.get<double>("elevation");

                std::string file = station.get<std::string>("file");
                auto f = cwd_dir / file;
                data.path = f.string();

                try
                {
                    auto filter_section = station.get_child("filter");

                    for (auto &jtr: filter_section)
                    {
                        auto filter_name = jtr.first.data();
                        auto cfg = jtr.second;

                        boost::shared_ptr<filter_base> filter = filter_factory::create(filter_name,cfg);
                        filter->init();

                        //save this filter to run later
                        //each station is associated with a list of filters to run
                        //metdata will own the filters after this
                        data.filters.push_back(filter);

                    }

                } catch (pt::ptree_bad_path &e)
                {
                    // no filters
                }

                ascii_data.push_back(data);

            }
        }
        _metdata->load_from_ascii(ascii_data, _global->_utc_offset);
        nstations = _metdata->nstations();
    }

    // based on config file the loaded metdata might be trimmed so we need to wait until
    // determine_starend_ts_forcing call to populate any global structs of star/end time.

    SPDLOG_DEBUG("Found # stations = {}", nstations);
    if(nstations == 0)
    {
        CHM_THROW_EXCEPTION(forcing_error, "No input forcing files found!");
    }

    SPDLOG_DEBUG("Finished reading stations. Took {} s", c.toc<s>());

}
void core::determine_startend_ts_forcing()
{

    if (_metdata->nstations() == 0)
    {
        CHM_THROW_EXCEPTION(forcing_error, "no stations");
    }

    auto start_time = _metdata->start_time();
    auto end_time = _metdata->end_time();

    //if the user gives both a custom start and end time, don't bother with this.
    //this does assume the user knows what they are doing and subsets a time period that is correct
    //if they do not, the sanity check after this period that double checks consistent timestepping will catch it
    if(!_start_ts && !_end_ts)
    {
        // find the latest start time and the earliest end time
        std::tie(start_time, end_time) = _metdata->start_end_time();
    }

    // If we ended on time T, restart from T+1. T+1 is written out to attr, so we can just start from this
    if(_checkpoint_opts.load_from_checkpoint)
    {
        size_t t = 0;
        _checkpoint_opts.in_savestate.get_ncfile().getAtt("restart_time_sec").getValues(&t);
        _start_ts = new boost::posix_time::ptime(boost::posix_time::from_time_t(t));

        SPDLOG_WARN("Loading from checkpoint. Overriding start time to match. New start time = {}", boost::posix_time::to_simple_string(*_start_ts));
    }

    if (!_start_ts)
    {
        _start_ts = new boost::posix_time::ptime(start_time);
    }
    else if(*_start_ts < start_time)
    {
        std::stringstream ss;
        ss << "User specified start time starts before the most continuous start time of all input forcing files . ";
        ss << " User start date: " << *_start_ts;
        ss << " Continous start date: " << start_time;
        CHM_THROW_EXCEPTION(model_init_error,ss.str());
    }
    if (!_end_ts)
    {
        _end_ts = new boost::posix_time::ptime(end_time);
    }
    else if(*_end_ts > end_time)
    {
        std::stringstream ss;
        ss << "User specified end time ends after the most continuous end time of all input forcing files . ";
        ss << " User end date: " << *_end_ts;
        ss << " Continous end date: " << end_time;
        CHM_THROW_EXCEPTION(model_init_error, ss.str());
    }


    if( *_start_ts > *_end_ts)
    {
        std::stringstream ss;
        ss << "Start time is after endtime. ";
        ss << " Start date: " << *_start_ts;
        ss << " End date: " << *_end_ts;
        CHM_THROW_EXCEPTION(model_init_error,ss.str());
    }


    //figure out what our timestepping is. Needs to happen before we subset as we may end up with only 1 timestep
    _global->_dt = _metdata->dt_seconds();
    SPDLOG_DEBUG("model dt = {}s", _global->dt());


    SPDLOG_DEBUG("Subsetting station timeseries");
    _metdata->subset(*_start_ts, *_end_ts);


    //ensure all the stations have the same start and end times
    // per-timestep agreement happens during runtime.
    _metdata->check_ts_consistency();

    // write this into global
    _global->_n_timestep = _metdata->n_timestep();
}
void core::config_parameters(pt::ptree &value)
{
    SPDLOG_DEBUG("Found parameter mapping section");

    //replace any references to external files with the file contents
    //we can't use _find_and_insert_subjson(value); as this needs to be handled slightly differently
    for(auto &itr : value)
    {
        if(itr.second.data().find(".json") != std::string::npos)
        {
            auto dir =  cwd_dir / itr.second.data();

            auto cfg = read_json(dir.string());
            std::string param_name = itr.first.data();
            //replace the string config name with that config file
            value.put_child(param_name, cfg);
        }
    }

    this->_global->parameters = value.get_child(""); //get root

}

bool core::config_meshes( pt::ptree &value)
{
    SPDLOG_DEBUG("Found meshes sections");

    _mesh = boost::make_shared<triangulation>();

    _mesh->_global = _global;

    _find_and_insert_subjson(value);

    _mesh_path = value.get<std::string>("mesh");
    SPDLOG_DEBUG("Found mesh:{}",_mesh_path);
    _mesh_path = (cwd_dir / _mesh_path).string();

    auto mesh_file_extension = boost::filesystem::path(_mesh_path).extension().string();

    bool is_partition = false;
    // only need to look for param + ic if we aren't loading a partition
    std::vector<std::string> param_file_paths;
    std::vector<std::string> initial_condition_file_paths;
    if(mesh_file_extension != ".partition")
    {

        // Paths for parameter files
        try
        {
            for(auto &itr : value.get_child("parameters"))
            {
                auto param_file = itr.second.data();
                SPDLOG_DEBUG("Found parameter file: {}",param_file);
                param_file_paths.push_back( (cwd_dir / itr.second.data()).string() );
            }
        }
        catch(pt::ptree_bad_path &e)
        {
            SPDLOG_DEBUG("No addtional parameters found in mesh section.");
        }

        // Paths for initial condition files
        try
        {
            for(auto &itr : value.get_child("initial_conditions"))
            {
                auto ic_file = itr.second.data();
                SPDLOG_DEBUG("Found initial condition file: {}",ic_file);
                initial_condition_file_paths.push_back( (cwd_dir / itr.second.data()).string() );
            }
        }
        catch(pt::ptree_bad_path &e)
        {
            SPDLOG_DEBUG("No addtional initial conditions found in mesh section.");
        }
    }


    // Ensure all files are HDF5 in multiprocess MPI runs
#ifdef USE_MPI
    if ( _comm_world.size() > 1 )
    {
        bool h5_or_part = mesh_file_extension != ".h5" || mesh_file_extension != ".partition";
        if(!h5_or_part)
        {
            CHM_THROW_EXCEPTION(mesh_error, "MPI multiprocess run requires hdf5 mesh.\n\n    Run the serial hdf5 conversion tool\n\n");
        }

        // only check the params and ics if we aren't using a partition file
        if(mesh_file_extension != ".partition")
        {

            for(const auto& it : param_file_paths)
            {
                auto extension = boost::filesystem::path(it).extension();
                if(extension != ".h5")
                {
                    CHM_THROW_EXCEPTION(mesh_error, "MPI multiprocess run requires hdf5 parameter files: " + it + "\n\n    Run the serial hdf5 conversion tool\n\n");
                }
            }
            for(const auto& it : initial_condition_file_paths)
            {
                auto extension = boost::filesystem::path(it).extension();
                if(extension != ".h5")
                {
                    CHM_THROW_EXCEPTION(mesh_error, "MPI multiprocess run requires hdf5 initial condition files: " + it + "\n\n    Run the serial hdf5 conversion tool\n\n");
                }
            }
        }

    }
#endif

    //we need to let the mesh know about any parameters the modules will provide so they can be correctly build into the static hashmaps
    for(auto& p : _provided_parameters) {
        _mesh->_parameters.insert(p);
    }

    _provided_parameters = _mesh->parameters();

    // Before we read the mesh, we need to know if we are geographic or projected so we can hook up all the distance functions
    // correctly. So for either json or hdf5 we check, hook up the functions, then proceed to the main load which can assume the functions are available
    bool is_geographic = check_is_geographic(_mesh_path);

    _global->_is_geographic = is_geographic; // save it here so modules can determine if this is true
    if( is_geographic)
    {
        math::gis::point_from_bearing = & math::gis::point_from_bearing_latlong;
        math::gis::distance = &math::gis::distance_latlong;
    }
    else
    {
        math::gis::point_from_bearing = &math::gis::point_from_bearing_UTM;
        math::gis::distance = &math::gis::distance_UTM;
    }

    ////////////////////////////////////////////////////////////
    // Actually read the mesh, parameter and ic data here
    ////////////////////////////////////////////////////////////
    if(mesh_file_extension == ".h5")
    {
        _mesh->from_hdf5(_mesh_path, param_file_paths, initial_condition_file_paths);
    }
    else if(mesh_file_extension == ".partition")
    {
#ifndef USE_MPI
        CHM_THROW_EXCEPTION(config_error, "Using a partitioned mesh requires enabling MPI during CHM configure and build.");
#endif
        is_partition = true;
        _mesh->from_partitioned_hdf5(_mesh_path);



    }
    else  // Assume anything that is NOT h5 is json
    {
        auto mesh = read_json(_mesh_path);

        //mesh will have been loaded by the geographic check so don't re load it here
      bool triarea_found = false;

      // Parameter files
      for (const auto& param_file : param_file_paths)
      {
	    pt::ptree param_json = read_json(param_file);

            for(auto& ktr : param_json)
            {
                //use put to ensure there are no duplicate parameters...
                std::string key = ktr.first.data();
                mesh.put_child( "parameters." + key ,ktr.second);

                if( key == "area")
                    triarea_found = true;

                SPDLOG_DEBUG("Inserted parameter {} into the config tree.", ktr.first.data());
	    }
      }

        if(is_geographic && !triarea_found)
        {
            CHM_THROW_EXCEPTION(mesh_error, "Geographic meshes require the triangle area be present in a .param file. Please include this.");
        }

      // Initial condition files
      for(const auto& ic_file : initial_condition_file_paths)
      {
            pt::ptree ic_json = read_json(ic_file);

            for(auto& ktr : ic_json)
            {
                //use put to ensure there are no duplciate parameters...
                std::string key = ktr.first.data();
                mesh.put_child( "initial_conditions." + key ,ktr.second);
                SPDLOG_DEBUG("Inserted initial condition {} into the config tree.",  ktr.first.data());
            }
      }

      _mesh->from_json(mesh);

    }

    if (_mesh->size_local_faces() == 0)
    {
      CHM_THROW_EXCEPTION(mesh_error, "Mesh size = 0!");
    }


    return is_partition;
}

void core::config_output(pt::ptree &value)
{
    SPDLOG_DEBUG("Found output section");
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkStringArray> labels = vtkSmartPointer<vtkStringArray>::New();
    labels->SetName("Point output name");

    size_t ID = 0;
    auto pts_dir = "points";
    boost::filesystem::path pts_path;

    auto output_dir = value.get<std::string>("output_dir","output");
    output_folder_path = cwd_dir / output_dir;
    boost::filesystem::create_directories(output_folder_path);

    // Create empty folders /points/
    pts_path = output_folder_path / pts_dir;
    boost::filesystem::create_directories(pts_path);

    _find_and_insert_subjson(value);

    for (auto &itr : value)
    {
        output_info out;
        std::string out_type = itr.first.data();
        if (out_type == "output_dir") //skip this
        {
            continue;
        }

        if (out_type == "mesh")
        {
            CHM_THROW_EXCEPTION(config_error,
                "output.mesh is now removed in favour of output.vtu or output.ugrid\n"
                "Please see\n\thttps://chm.readthedocs.io/en/develop/configuration.html#output \n for more information.");
        }

        if ((out_type !="vtu") &&
            (out_type !="ugrid"))  // anything else *should* be a time series*......
        {
            out.type = output_info::time_series;
            out.name = out_type; //station name = key name

            std::string fname;
            try
            {
                fname = itr.second.get<std::string>("file");
            }
            catch(const pt::ptree_error &e)
            {
                CHM_THROW_EXCEPTION(forcing_error,"Missing output filename for " + out.name);
            }
            auto f = pts_path / fname;
            out.fname = f.string();

            try
            {
                out.longitude = itr.second.get<double>("longitude");
                out.latitude = itr.second.get<double>("latitude");
            }
            catch(const pt::ptree_error &e)
            {
                CHM_THROW_EXCEPTION(forcing_error, "Output point " + out.name + " is missing latitude and/or longitude.");
            }


            if( (out.latitude > 90 || out.latitude < -90) ||
                (out.longitude > 180 || out.longitude < -180) )
            {
                CHM_THROW_EXCEPTION(forcing_error, "Output " + out.name + " coordinate is invalid.");
            }

            //project mesh, need to convert the input lat/long into the coordinate system our mesh is in
            if(!_mesh->is_geographic())
            {

                OGRSpatialReference insrs;
                insrs.SetWellKnownGeogCS("WGS84");
                insrs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER); //enforce ingoing as x y

                OGRSpatialReference outsrs;
                outsrs.importFromProj4(_mesh->proj4().c_str());
                outsrs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER); //enforce outgoing as x y

                out.x = out.longitude;
                out.y = out.latitude;

                OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&insrs, &outsrs);

                //CRS created with the “EPSG:4326” or “WGS84” strings use the latitude first, longitude second axis order.
                if(!coordTrans->Transform(1, &out.x, &out.y))
                {
                    CHM_THROW_EXCEPTION(forcing_error, "Output=" + out.name + ": unable to convert coordinates to mesh format.");
                }

                OGRCoordinateTransformation::DestroyCT(coordTrans);
            }

            if(!_mesh->is_geographic())
                out.face = _mesh->locate_face(out.x, out.y);
            else
                out.face = _mesh->locate_face(out.longitude, out.latitude);

            if(out.face != nullptr)
            {
                // In MPI mode, we might not have the output triangle on this node, so we just have to fail gracefully and hope that the other nodes have this.
                // in the future we need a comms here to check if the nodes successfully figured out the output points
                // for now, if we don't have it, print an error but keep going

                //set the point to be the center of the triangle that the output point lies on
                points->InsertNextPoint(out.face->get_x(), out.face->get_y(), out.face->get_z());
                labels->InsertNextValue(out.name);

                SPDLOG_DEBUG("Triangle geometry for output triangle = {} slope: {} aspect:{}",
                              out_type, out.face->slope() * 180. / M_PI, out.face->aspect() * 180. / M_PI);

                out.face->_debug_name = out.name; //out_type holds the station name
                out.face->_debug_ID = ID;
                ++ID;
            }
            else
            {
                SPDLOG_WARN("Output point {} was not found in this rank's mesh", out.name);
            }
        }
        else if (out_type == "vtu" || out_type == "ugrid")
        {
            out.type = output_info::mesh;

            out.name = out_type;
            out.base_name = itr.second.get<std::string>("base_name", output_dir);

            auto write_all_parameters = itr.second.get_optional<bool>("write_all_parameters");
            if (!write_all_parameters)
            {
                write_all_parameters = itr.second.get_optional<bool>("write_parameters");
                if (write_all_parameters)
                {
                    SPDLOG_WARN("Output {} uses deprecated write_parameters; use write_all_parameters instead.", out.name);
                }
            }
            _mesh->write_param_to_output(write_all_parameters.value_or(true));

            try
            {
                for (auto &jtr: itr.second.get_child("variables"))
                {
                    out.variables.insert(jtr.second.data());
                }
            }
            catch (pt::ptree_bad_path &e)
            {
                SPDLOG_WARN("Writing all variables to output mesh");
            }

            try
            {
                for (auto &jtr: itr.second.get_child("output_parameters"))
                {
                    out.output_parameters.insert(jtr.second.data());
                }
            }
            catch (pt::ptree_bad_path &e)
            {
                // defaults handled in triangulation
            }

            if (out_type == "vtu")
            {
                auto msh_path = output_folder_path / "vtu";
                boost::filesystem::create_directories(msh_path);

                auto f = msh_path / out.base_name;
                boost::filesystem::create_directories(f.parent_path());
                out.fname = f.string();
                out.mesh_output_formats = output_info::mesh_outputs::vtu;
            }
            else
            {
                auto format = itr.second.get<std::string>("format", "netcdf");
                bool use_zarr = false;
                std::string extension = ".nc";

                if (format == "zarr")
                {
                    use_zarr = true;
                    extension = ".zarr";
                }
                else if (format != "netcdf" && format != "nc")
                {
                    CHM_THROW_EXCEPTION(config_error, "Unknown ugrid format: " + format);
                }

                boost::filesystem::path ugrid_path = output_folder_path / (out.base_name + extension);
                out.mesh_output_formats = output_info::mesh_outputs::ugrid;
                out.fname = ugrid_path.string();

                out.writer = boost::make_shared<ugrid_writer>(
                    _mesh,
                    _global,
                    _mesh->write_param_to_output(),
                    out.fname,
                    use_zarr);

                boost::get<boost::shared_ptr<ugrid_writer>>(out.writer)->compress = itr.second.get<bool>("compress", true);
                boost::get<boost::shared_ptr<ugrid_writer>>(out.writer)->bitgroom = itr.second.get<bool>("bitgroom", true);

            }


	    // Set option for writing ghost neighbor data, defaults to not
            // only possible for vtu outputs
            auto write_ghost = itr.second.get("write_ghost_neighbors",false);
            if (write_ghost && out.mesh_output_formats == output_info::mesh_outputs::ugrid )
            {
                CHM_THROW_EXCEPTION(config_error, "ugrid output cannot have write_ghost_neighbors=true");
            }
            _mesh->write_ghost_neighbors_to_vtu(write_ghost) ;

            out.frequency = itr.second.get_optional<size_t>("frequency"); //defaults to every timestep
            out.rotate_frequency = itr.second.get_optional<size_t>("rotate_frequency"); //defaults to never

            if (out.rotate_frequency)
            {
                SPDLOG_DEBUG("Creating new UGRID output every {} timestep", *out.rotate_frequency);
            }

            out.only_last_n = itr.second.get_optional<size_t>("only_last_n");

            if(out.frequency && out.only_last_n)
            {
                SPDLOG_WARN("Only only_last_n output option will be used");
                out.frequency.reset();
            }

            if (out.mesh_output_formats == output_info::mesh_outputs::ugrid)
            {
                boost::get<boost::shared_ptr<ugrid_writer>>(out.writer)
                    ->set_output_cadence(out.frequency, out.only_last_n);
            }

            auto specific_datetime = itr.second.get_optional<std::string>("specific_datetime");
            if(specific_datetime)
            {
                out.specific_datetime = boost::posix_time::ptime(boost::posix_time::from_iso_string(*specific_datetime));
            }

            auto specific_time = itr.second.get_optional<std::string>("specific_time");
            if(specific_time)
            {
                // ptime needs a real date, but we will only ever be checking the minute and hour
                auto time = "3000-01-01 " + *specific_time + ":00";
                out.specific_time = boost::posix_time::ptime(boost::posix_time::time_from_string(time));
            }


            out.list_outputs();

        } else
        {
            SPDLOG_WARN("Unknown output type: {}",itr.second.data());
        }

        //ensure we aren't adding timeseries outputs with invalid faces bc of MPI
        if(  out.type == output_info::time_series &&
             out.face == nullptr)
        {
            // we will pass for now, but ultimiately we should have done an MPI comms and
            // check if we are missing an output

            SPDLOG_WARN("There is currently no check if all the MPI ranks correctly find the output timeseries triangle. "
                           "If you are missing output, ensure that all the output points are within the domain.");


        }else
        {
            _outputs.push_back(out);
        }
    }

    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(labels);
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    //output this to the same folder as the points are written out to
    std::string rank = "."+ std::to_string(_comm_world.rank());

    auto f = pts_path / ("output_points"+rank+".vtp");
    writer->SetFileName(f.string().c_str());
    #if VTK_MAJOR_VERSION <= 5
        writer->SetInput(polydata);
    #else
        writer->SetInputData(polydata);
    #endif

    writer->Write();

}

void core::config_global( pt::ptree &value)
{
    CHM_THROW_EXCEPTION(no_modules_defined,"Global section is now deprecated. This does nothing, please remove from config file.");
}

core::cmdl_opt core::config_cmdl_options(int argc, char **argv)
{


    std::string config_file = "";
    std::string start;
    std::string end;

    bool legacy_log=false;

    po::options_description desc("Allowed options.");
    desc.add_options()
            ("help", "This message")
            ("version,v", "Version number")
            ("legacy-log,l",po::value<bool>(&legacy_log), "Legacy log output writes to CHM.log instead of a timestamped CHM_*.log in log/")
            ("config-file,f", po::value<std::string>(&config_file),
             "Configuration file to use. Can be passed without --config-file [-f] as well ")
            ("config,c", po::value<std::vector<std::string>>(), "Specifies a configuration parameter."
                    "This can over-ride existing values."
                    "The value is specified with a fully qualified config path. Does not support list values []"
                    "For example:\n"
                    "-c config.Harder_precip_phase.const.b:1.5 -c config.debug.debug_level:\"error\"")
            ("remove,r", po::value<std::vector<std::string>>(), "Removes a configuration paramter."
                    " Removals are processed after any --config paramters are parsed, so -r will override -c. "
                    "Does not support remove of individual list [] itmes. For example:"
                    "-c nproc:2 -r nproc\n "
                    "will result in nproc being removed.")
            ("remove-module,d", po::value<std::vector<std::string>>(), "Removes a module."
                    " Removals are processed after any --config paramters are parsed, so -d will override -c. ")
            ("add-module,m", po::value<std::vector<std::string>>(), "Adds a module.");



    //allow for specifgying config w/o --config
    po::positional_options_description p;
    p.add("config-file", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        cout << desc << std::endl;
        CHM_THROW_EXCEPTION(chm_done,"done");
    }
    else if (vm.count("version"))
    {
        cout << version << std::endl;
        CHM_THROW_EXCEPTION(chm_done,"done");
    }

    std::vector<std::pair<std::string, std::string>> config_extra;
    if (vm.count("config"))
    {
        boost::char_separator<char> sep(":");
        for (auto &itr : vm["config"].as<std::vector<std::string>>())
        {
            boost::tokenizer<boost::char_separator<char>> tok(itr, sep);
            std::vector<std::string> v;

            for (auto &jtr : tok)
                v.push_back(jtr);

            if (v.size() != 2)
            {
                CHM_THROW_EXCEPTION(io_error, "Config value of " + itr + " is invalid.");
            }

            std::pair<std::string, std::string> override;
            override.first = v[0];
            override.second = v[1];
            config_extra.push_back(override);
        }
    }

    std::vector<std::string> rm_config_extra;
    if (vm.count("remove"))
        rm_config_extra = vm["remove"].as<std::vector<std::string>>();

    std::vector<std::string> remove_module;
    if (vm.count("remove-module"))
        remove_module = vm["remove-module"].as<std::vector<std::string>>();

    std::vector<std::string> add_module;
    if (vm.count("add-module"))
        add_module = vm["add-module"].as<std::vector<std::string>>();

    if (!vm.count("config-file"))
    {
        SPDLOG_ERROR("Configuration file required.");
        cout << desc << std::endl;
        exit(1);
    }


    return boost::make_tuple(config_file, //0
                             config_extra,  //1
                             rm_config_extra,  //2
                             remove_module,  //3
                             add_module, //4
                             legacy_log); //5
}

void core::init(int argc, char **argv)
{
    //default logging level
    _log_level = debug;

    auto log_start_time = boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());

    //get any command line options
    auto cmdl_options = config_cmdl_options(argc, argv);

    //see if we are getting passed a path for our config file. If so, we need to append this path to
    //any config files we are about to read in the main config file.

    //do this here though as we need cwd_dir for the log output path
    boost::filesystem::path path(cmdl_options.get<0>());
    cwd_dir = boost::filesystem::current_path();

    std::string log_dir = "log";
    auto log_path = cwd_dir / log_dir;

    //output a unique logfile for each mpi rank
    std::string rank = "";
    rank = "."+std::to_string(_comm_world.rank());

    std::string log_name = "CHM_" + log_start_time + rank + ".log";

    if (_comm_world.rank() == 0)
        boost::filesystem::create_directories(log_path);
    _comm_world.barrier();


    log_file_path = log_path / log_name;
    log_name =  log_file_path.string();

    auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    console_sink->set_level(spdlog::level::debug);
    console_sink->set_pattern("[%n] [%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");

    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_name, true);
    file_sink->set_pattern("[%n] [%Y-%m-%d %H:%M:%S.%e] [%s:%!:%#] [%^%l%$] %v");;
    file_sink->set_level(spdlog::level::debug);

    spdlog::sinks_init_list sink_list = { file_sink, console_sink };

    // make default logger
    // logger name is the current MPI rank
    spdlog::set_default_logger(std::make_shared<spdlog::logger>("rank " + std::to_string(_comm_world.rank()), spdlog::sinks_init_list({console_sink, file_sink})));
//    spdlog::set_pattern("[%P] [%s:%!:%#] [%^%l%$] %v");
    spdlog::set_level(spdlog::level::debug); // this is needed as well as the above levels
    spdlog::flush_on(spdlog::level::debug);

    SPDLOG_DEBUG("Logger initialized. Writing to cout and {}",  log_name);


    _global = boost::make_shared<global>();

    // This needs to be set so that underflows in gsl math
    // computations are not treated as errors. i.e.,
    //    gsl: expint.c:363: ERROR: underflow
    //    Default GSL error handler invoked.
    // is not what we want.
    gsl_set_error_handler_off();


    SPDLOG_DEBUG( "Built with MPI support,    #processes = {}", _comm_world.size());


#ifdef _OPENMP
    SPDLOG_DEBUG( "Built with OpenMP support, #threads   = {}", omp_get_max_threads());
#endif

    SPDLOG_DEBUG("PID={}",getpid());

    _hpc_scheduler_info.detect();

    pt::ptree cfg;
    try
    {
        SPDLOG_DEBUG("Reading configuration file {}", cmdl_options.get<0>());
        cfg = read_json(cmdl_options.get<0>());

        /*
         * The module config section is optional, but if it exists, we need it to
         * setup the modules, so check if it exists
         */
        //for each module: config pair
        for (auto &itr: cfg.get_child("config"))
        {
            pt::ptree module_config;

            //should we try to load a config file? if second.data() is a substree it's string is "" so this will properly work
            if(itr.second.data().find(".json") != std::string::npos)
            {
                //load the module config into it's own ptree
                auto dir =  cwd_dir / itr.second.data();

                module_config = read_json(dir.string());

                std::string module_name = itr.first.data();
                //replace the string config name with that config file
                cfg.put_child("config." + module_name, module_config);
            }

        }

    }
    catch (pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG( "Optional section Module config not found" );
    }
    catch (pt::json_parser_error &e)
    {
        CHM_THROW_EXCEPTION(config_error, "Error reading file: " + cmdl_options.get<0>() + " on line: " + std::to_string(e.line()) + " with error: " +
                e.message());
    }



    //now apply any override or extra configuration parameters from the command line
    for (auto &itr : cmdl_options.get<1>())
    {
        //check if we are overwriting something
        try
        {
            auto value = cfg.get<std::string>(itr.first);
            SPDLOG_WARN( "Overwriting {} = {} with {} = {}  ", itr.first, value, itr.first, itr.second );

        } catch (pt::ptree_bad_path &e)
        {
            SPDLOG_DEBUG("Inserting new config {} = {}", itr.first, itr.second);
        }

        cfg.put(itr.first, itr.second);
    }



    //now apply any removal parameters from the command line
    for (auto &itr:cmdl_options.get<2>())
    {
        //check if we are overwriting something
        try
        {
            auto value = cfg.get<std::string>(itr);

            std::size_t found = itr.rfind(".");

            //if we have a subkey
            if (found != std::string::npos)
            {
                std::string parent = itr.substr(0, found);
                std::string child = itr.substr(found + 1, itr.length() - found + 1);

                cfg.get_child(parent).erase(child);
            }
            else //otherwise just blow it away.
            {
                cfg.erase(itr);
            }


            SPDLOG_DEBUG( "Removing {}", itr);

        } catch (pt::ptree_bad_path &e)
        {
            SPDLOG_DEBUG("No value {} to remove", itr);
        }

    }

    //remove and add modules
    /*
     * We expect the following sections:
     *  modules <- These don't get the mesh in the ctor, everything has to be in init, so doesn't mess with the MPI calls
     *  meshes  <- Will redefine the elements to be on a per-node (MPI mode) basis
     *  forcing
     * The rest may be optional, and will override the defaults.
     */
    config_modules(cfg.get_child("modules"), cfg.get_child("config"), cmdl_options.get<3>(), cmdl_options.get<4>());

    //update the internal config ptree so when we right it out for auditing
    cfg.erase( "modules") ;

    pt::ptree array;
    for(auto& itr : get_active_module_list())
        array.push_back(pt::ptree::value_type("", itr.first->ID));

    cfg.put_child("modules", array);

    // This has the delayed param load enabled, so mesh path is saved to _mesh_path which is used to load
    // the params latter
     bool ispart = config_meshes(cfg.get_child("meshes")); // this must come before forcing, as meshes initializes the required distance functions based on geographic/utm meshes


    //output should come before forcing, controls if we should output the vtp file of station locations
    try
    {
        config_output(cfg.get_child("output"));
    } catch (pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG( "Optional section Output not found");
    }


    // needs both the mesh loaded and output folder setup
    boost::filesystem::create_directories(output_folder_path / "mesh_boundingbox");
    auto f = output_folder_path / "mesh_boundingbox" / fmt::format("mesh_bbox_{}.geojson", _comm_world.rank());
    _mesh->write_bbox_geojson(f.string());

    config_forcing(cfg.get_child("forcing"));

    /*
     * We can expect the following sections to be optional.
     */
    try
    {
        config_options(cfg.get_child("option"));
    } catch (pt::ptree_bad_path &e)
    {
        CHM_THROW_EXCEPTION(config_error, "The configuration option section is missing and is required.");
    }


    // Now the forcing and mesh are loaded, assign each face the station lists
    populate_face_station_lists();

    //load the parameters now that the station list has been pruned and we have a partitioned mesh
    if( ispart)
        _mesh->from_partitioned_hdf5(_mesh_path, true);

    boost::filesystem::path full_path(boost::filesystem::current_path());

    try
    {
        config_global(cfg.get_child("global"));
    } catch (pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG("Optional section Global not found");
    }

    try
    {
        config_module_overrides(cfg.get_child("remove_depency"));
    } catch (pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG( "Optional section remove_depency not found");
    }

    try
    {
        config_parameters(cfg.get_child("parameter_mapping"));
    } catch (pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG("Optional section parameter mapping not found");
    }

    try
    {
        config_checkpoint(cfg.get_child("checkpoint"));
    } catch (pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG("Optional section checkpoint not found");
    }

    //if we are to run in point mode, we need to remove all the input and outputs that aren't associated with point mode
    // we do this so the user doesn't have to modify a lot of the config file when going between point and dist mode.

    if(point_mode.enable)
    {
        SPDLOG_DEBUG("Running in point mode");
        //Each face knows which stations are closest to it and what it should use

        if(point_mode.use_specific_station && _metdata->is_netcdf())
        {
            CHM_THROW_EXCEPTION(model_init_error, "If a specific station is requested, this must be done using an ascii forcing file definition." );
        }
        if(point_mode.use_specific_station && _outputs.size() !=1)
        {
            CHM_THROW_EXCEPTION(model_init_error, "If a specific station is requested, there must be only 1 output location. This allows for a 1:1 mapping of station to output point." );
        }

        if(point_mode.use_specific_station)
        {
            std::unordered_set<std::string> remove_set;
            // remove everything but the one forcing that the user asked for
            for (auto& itr : _metdata->stations())
            {
                if (itr->ID() != point_mode.forcing)
                    remove_set.insert(itr->ID());
            }

            _metdata->prune_stations(remove_set);

            if ( _metdata->nstations() != 1 )
            {
                CHM_THROW_EXCEPTION(model_init_error, "The number of stations is " + std::to_string(_metdata->nstations()) + ". In point mode must be equal to exactly 1" );
            }
        }

        if(_metdata->is_netcdf())
        {
            SPDLOG_DEBUG("Using the following input nc grid cells as forcing:");
            for(auto& s:_outputs.at(0).face->stations())
            {
                SPDLOG_DEBUG( "\tx={} y={} lat={} lon={}", s->_nc_x, s->_nc_y, s->y(),s->x());
            }
        }
    }

    // output a full dump of the cfg, after all modifications, to the output directory
    pt::json_parser::write_json((output_folder_path / "config.json" ).string(), cfg);
    _cfg = cfg;

    SPDLOG_DEBUG("Finished initialization");

    SPDLOG_DEBUG("Determining module dependencies");
    _determine_module_dep();

    //now we know what outputs we have, and have ensure that's valid, we need to ensure the user hasn't asked to output
    // a variable that won't be created, otherwise this will segfault.

    for(auto& o : _outputs)
    {
        if( o.type != output_info::mesh )
            continue;

        for (auto& var : o.variables)
        {
            //check every output we requested against the global full list of provided outputs, bail if we don't find it
            if(!boost::algorithm::any_of_equal(_provided_var_module,var))
            {
                CHM_THROW_EXCEPTION(config_error, "Requested output " + var + " is not provided by any module.");
            }
        }
    }

    determine_startend_ts_forcing();

    //set interpolation algorithm
    _global->interp_algorithm = _interpolation_method;


    //setup output timeseries sinks

    for (auto &itr : _outputs)
    {
        if (itr.type == output_info::output_type::time_series)
        {
            itr.ts.init(_provided_var_module, _metdata->start_time(), _metdata->end_time(), _metdata->dt());
        }
    }

    if(point_mode.enable)
    {
        for(const auto& itr:_chunked_modules)
        {
            for(auto jtr:itr)
            {
                if(jtr->parallel_type() == module_base::parallel::domain)
                {
                    CHM_THROW_EXCEPTION(model_init_error, "Domain parallel module are being run in point-mode.");
                }
            }
        }
    }

    SPDLOG_DEBUG( "Allocating face variable storage");

    //we are going to make the assumption that every module can store face data.
    // However if this is onerous we can add a flag to the modules later
    std::set<std::string> module_list;
    for(auto& itr: _modules)
    {
        module_list.insert(itr.first->ID);
    }

    // only init on the faces we are going to compute on
    // seemed like a good idea but module init has no way to know that it shouldn't run on part of the mesh
    // todo: revisit only init on parts of the mesh
    if(point_mode.enable)
    {
        SPDLOG_DEBUG("Initialzing faces for point_mode only");
        std::vector<mesh_elem> faces_to_init;

        for (auto &itr : _outputs)
        {
            if (itr.type == output_info::output_type::time_series)
            {
                faces_to_init.push_back(itr.face);
            }
        }
        _mesh->prune_faces(faces_to_init);
        SPDLOG_DEBUG("Mesh now has #faces = {}",_mesh->size_local_faces());
    }

    // module provided params are initialized when the mesh loads as we also have to handle the mesh params at the same time
    _mesh->init_face_data(_provided_var_module, _provided_var_vector, module_list);

    timer c;
    SPDLOG_DEBUG("Running init() for each module");
    c.tic();


    for (auto& itr : _modules)
    {
      spdlog::debug("\t{}", itr.first->ID);
      itr.first->init(_mesh);
    }



    SPDLOG_DEBUG("Took {}ms", c.toc<ms>());

    //we do this here now because init is allowing a module to chance its mind and declare itself
    // data parallel or domain parallel after the fact.
    _schedule_modules();

//load a checkpoint as the last thing we do before a run
    if(_checkpoint_opts.load_from_checkpoint  )
    {
        _global->_from_checkpoint = true;
        SPDLOG_DEBUG("Loading from checkpoint");
        c.tic();
        for (auto &itr : _chunked_modules)
        {
            //module calls
            for (auto &jtr : itr)
            {
                jtr->load_checkpoint(_mesh, _checkpoint_opts.in_savestate);
            }
        }

        SPDLOG_DEBUG("Done loading snapshot [ {}s ]", c.toc<s>());
    }
}

void core::_find_and_insert_subjson(pt::ptree& value)
{
    std::vector<std::string> keys_to_remove; //probably can't remove a node without invalidating the iterator. so mark it and remove it later.
    //replace any references to external files with the file contents
    for(auto &itr : value)
    {
        if(itr.second.data().find(".json") != std::string::npos)
        {
            auto dir =  cwd_dir / itr.second.data();

            auto cfg = read_json(dir.string());
            std::string key = itr.first.data();
            keys_to_remove.push_back(key);

            for (auto &jtr : cfg)
            {
                //If an external file is provided, don't allow further sub json files. That is a rabbit hole that isn't worth dealing with.
                if (jtr.second.data().find(".json") != std::string::npos)
                {
                    CHM_THROW_EXCEPTION(forcing_error, "An included json file cannot contain json sub-files.");
                }

                std::string point_name = jtr.first.data();
                //replace the string config name with that config file
                value.put_child(point_name, jtr.second);
            }
        }
    }

    for(auto& itr : keys_to_remove)
    {
        value.erase(itr);
    }
}

void core::_determine_module_dep()
{

    size_t size = _modules.size();

    //init a graph of the required size == num of modules
    Graph g(size);


    std::set<std::string> graphviz_vars;

    bool overwrite_found = false;
    //make sure modules cannot override another module's output
    for (auto &m1_pair : _modules)
    {
      auto m1 = m1_pair.first;
      // Get vector of module m1's variable names
      auto m1_provides_var_names = m1->get_variable_names_from_collection(*(m1->provides()));;

        for (auto &m2_pair : _modules)
        {
	  auto m2 = m2_pair.first;
            if(m1->ID == m2->ID )
                continue;

      // Get vector of current module's variable names
      auto m2_provides_var_names = m2->get_variable_names_from_collection(*(m2->provides()));;;

           if( std::find(m1->conflicts()->begin(),
                   m1->conflicts()->end(),
                   m2->ID) != m1->conflicts()->end())
           {
               CHM_THROW_EXCEPTION(module_error,  "Module " + m1->ID + " explicitly conflicts with " + m2->ID);
           }


            auto itr = std::find_first_of (m1_provides_var_names.begin(), m1_provides_var_names.end(),
                                           m2_provides_var_names.begin(), m2_provides_var_names.end());
            if( itr != m1_provides_var_names.end() )
            {
                overwrite_found=true;
                SPDLOG_ERROR("Module {} and {} both provide variable {} ", m1->ID, m2->ID, *itr);
            }
        }
    }

    if(overwrite_found)
    {
        CHM_THROW_EXCEPTION(module_error,
                            "A module's provides overwrites another module's provides. This is not allowed.");
    }

    //build a list of variables provided by the met files, culling duplicate variables from multiple stations.

    auto vars = _metdata->list_variables();
    _provided_var_met_files.insert(vars.begin(), vars.end());


    //loop through each module
    for (auto &module_pair : _modules)
    {
      auto module = module_pair.first;

      // Get vector of current module's variable names
      auto mod_provides_var_names = module->get_variable_names_from_collection(*(module->provides()));

        //Generate a  list of all variables,, provided from this module, and append to the total list, culling duplicates between modules.
        _provided_var_module.insert(mod_provides_var_names.begin(), mod_provides_var_names.end());

//        vertex v;
//        v.name = module.first->ID;
//        boost::add_vertex(v,g);
        g[module->IDnum].name = module->ID;
//        output_graph << module.first->IDnum << " [label=\"" << module.first->ID  << "\"];" << std::endl;

        //names.at(module.first->IDnum) =  module.first->ID;

        //check intermodule depends
        if (module->depends()->size() == 0)  //check if this module requires dependencies
        {
            SPDLOG_DEBUG("Module [{}], no inter-module dependenices", module->ID);
        } else
        {
            SPDLOG_DEBUG("Module [{}], checking against...", module->ID);
        }


        //make a copy of the depends for this module. we will then count references to each dependency
        //this will allow us to know what is missing
        std::map<std::string, size_t> curr_mod_depends;
	auto curr_mod_depends_var_names = module->get_variable_names_from_collection(*(module->depends()));

        //populate a list of all the dependencies
        for (auto &itr : curr_mod_depends_var_names)
        {
            curr_mod_depends[itr] = 0; //ref count init to 0
        }

        //iterate over all the modules,
        for (auto &itr_pair : _modules)
        {
	  auto itr_module = itr_pair.first;

	  // Get vector of current module's variable names
	  auto itr_mod_provides_var_names = itr_module->get_variable_names_from_collection(*(itr_module->provides()));
            //don't check against our module
            if (module->ID.compare(itr_module->ID) != 0)
            {
                //loop through each required variable of our current module
                for (auto &depend_var : *(module->depends()))
                {
                    //LOG_DEBUG << "\t\t[" << itr_module->ID << "] looking for var=" << depend_var;

                    auto i = std::find(itr_mod_provides_var_names.begin(), itr_mod_provides_var_names.end(), depend_var.name);
                    if (i != itr_mod_provides_var_names.end()) //itr provides the variable we are looking for
                    {
                        SPDLOG_DEBUG("\t\tAdding edge between {} [{}] -> {} [{}] for var={}", module->ID, module->IDnum, itr_module->ID,itr_module->IDnum, *i);


                        //add the dependency from module -> itr, such that itr will come before module
                        edge e;
                        e.variable = *i;

                        //check if we need to ignore this edge as a result of user specific overrdige
                        //this will avoid a cycle is this exists.
                        bool ignore = false;
                        for (auto o : _overrides)
                        {
                            if (o.first == module->ID &&
                                o.second == itr_module->ID)
                            {
                                ignore = true;
                                SPDLOG_DEBUG("Skipped adding edge betweenn {} and {} for var={}  because of user override", o.first, o.second,*i);

                            }
                        }

                        if (!ignore)
                            boost::add_edge(itr_module->IDnum, module->IDnum, e, g);

                        //even if we ignore, inc our depencies so we don't fail later

                        //output_graph << itr_module->IDnum << "->" << module->IDnum << " [label=\"" << *i << "\"];" << std::endl;
                        curr_mod_depends[*i]++; //ref count our variable

                        //only link in the output if we aren't ignoring
                        if (!ignore)
                            graphviz_vars.insert(*i);
                    }
                }

                //loop through each required variable of our current module
                for (auto &optional_var : *(module->optionals()))
                {
                    //LOG_DEBUG << "\t\t[" << itr_module->ID << "] looking for var=" << depend_var;

		  auto itr_mod_provides_var_names = itr_module->get_variable_names_from_collection(*(itr_module->provides()));

                    auto i = std::find(itr_mod_provides_var_names.begin(), itr_mod_provides_var_names.end(), optional_var);
                    if (i != itr_mod_provides_var_names.end()) //itr provides the variable we are looking for
                    {
                        SPDLOG_DEBUG("\t\tAdding edge between {} [{}] -> {} [{}] for var={}", module->ID, module->IDnum, itr_module->ID,itr_module->IDnum, *i);


                        //add the dependency from module -> itr, such that itr will come before module
                        edge e;
                        e.variable = *i;

                        //check if we need to ignore this edge as a result of user specific override
                        //this will avoid a cycle if this exists.
                        bool ignore = false;
                        for (auto o : _overrides)
                        {
                            if (o.first == module->ID &&
                                o.second == itr_module->ID)
                            {
                                ignore = true;
                                SPDLOG_DEBUG("Skipped adding edge betweenn {} and {} for var={}  because of user override", o.first, o.second,*i);
                            }
                        }

                        if (!ignore)
                            boost::add_edge(itr_module->IDnum, module->IDnum, e, g);

                        //output_graph << itr_module->IDnum << "->" << module->IDnum << " [label=\"" << *i << "\"];" << std::endl;
                        //curr_mod_depends[*i]++; //ref count our variable

                        graphviz_vars.insert(*i);

                        module->set_optional_found(*i);
                    }
                }
            }
        }


        bool missing_depends = false;

        std::stringstream ss;
        //build a list of ALL missing variables before dying
        for (auto &itr : curr_mod_depends)
        {
            if (itr.second == 0)
            {
                ss << "Missing inter-module dependencies for module [" << module->ID << "]: " << itr.first << "\n";
                missing_depends = true;
            }

        }
        if (missing_depends)
        {
            SPDLOG_ERROR((ss.str()));
            CHM_THROW_EXCEPTION(module_error,   ss.str());
        }

        //check if our module has any met file dependenices
        if (module->depends_from_met()->size() == 0)
        {
            SPDLOG_DEBUG("Module [{}], no met file dependencies", module->ID);
        } else
        {
            SPDLOG_DEBUG("Module [{}] has met file dependencies",module->ID);
        }




        //check this modules met dependencies, bail if we are missing any.
        for (auto &depend_met_var : *(module->depends_from_met()))
        {
            auto i = std::find(_provided_var_met_files.begin(), _provided_var_met_files.end(), depend_met_var);
            if (i == _provided_var_met_files.end())
            {
                SPDLOG_ERROR("\t\t[]...[missing]",depend_met_var);
                CHM_THROW_EXCEPTION(module_error,  "Missing dependency for " + depend_met_var);
            }
            SPDLOG_DEBUG("\t\t{}...[ok]", depend_met_var);
        }

    }

    //great filter file for gvpr
    std::ofstream gvpr("filter.gvpr");;

    std::string font = "Helvetica";
    int fontsize = 11;

    std::string edge_str = "E[edgetype == \"%s\"] {\n color=\"/paired12/%i\";\n fontsize=%i;\n     fontname=\"%s\"\n }";
    int idx = 1;
    for (const auto& itr : graphviz_vars)
    {
        std::string edge = str_format(edge_str, itr.c_str(), idx, fontsize, font.c_str());
        idx++;
        gvpr << edge << std::endl;
    }

    gvpr.close();

    std::deque<int> topo_order;
    bool failed_toposort_call = false; // we need to still write the graph out if our topo sort fails
    try
    {
        boost::topological_sort(g, std::front_inserter(topo_order));
    }
    catch(...)
    {
        failed_toposort_call = true; // get ready to bail after write out the topo call
    }

    std::ostringstream ssdot;
    boost::write_graphviz(ssdot, g, boost::make_label_writer(boost::get(&vertex::name, g)),
                          make_edge_writer(boost::get(&edge::variable, g)));
    std::string dot(ssdot.str());
    size_t pos = dot.find("\n");

    if (pos == std::string::npos)
    {
        CHM_THROW_EXCEPTION(config_error, "Unable to generate dot file");
    }

    //skip past the newline
    pos++;

//    Insert the following to make the chart go right to left, landscape
//    rankdir=LR;
//    {
//        node [shape=plaintext, fontsize=16];
//        "Module execution order"->"";
//    }
    dot.insert(pos,
               "rankdir=LR;\n{\n\tnode [shape=plaintext, fontsize=16];\n\t\"Module execution order\"->\"\";\n}\nsplines=polyline;\n");

    std::ofstream file;
    file.open("modules.dot.tmp");
    file << dot;
    file.close();

    //http://stackoverflow.com/questions/8195642/graphviz-defining-more-defaults
    int ierr = std::system("gvpr -c -f filter.gvpr -o modules.dot modules.dot.tmp > /dev/null 2>&1"); CHK_SYSTEM_ERR(ierr);
    ierr = std::system("dot -Tpdf modules.dot -o modules.pdf > /dev/null 2>&1"); CHK_SYSTEM_ERR(ierr);
    ierr = std::remove("modules.dot.tmp"); CHK_SYSTEM_ERR(ierr);
    ierr = std::remove("filter.gvpr"); CHK_SYSTEM_ERR(ierr);

    // if we failed to do the topological sort above, now we can safety bail as we've written out the module.pdf
    if(failed_toposort_call)
    {
        CHM_THROW_EXCEPTION(config_error,  "Module graph must be a DAG. Please review modules.pdf to determine where the cycle occurred.");
    }

    std::stringstream ss;
    size_t order = 0;
    for (std::deque<int>::const_iterator i = topo_order.begin(); i != topo_order.end(); ++i)
    {
        ss << _modules.at(*i).first->ID << "->";
        _modules.at(*i).second = order;
        order++;
    }

    std::string s = ss.str();
    SPDLOG_DEBUG("Build order: {}",s.substr(0, s.length() - 2));


    //sort ascending based on make order number
    std::sort(_modules.begin(), _modules.end(),
              [](const std::pair<module, size_t> &a, const std::pair<module, size_t> &b) -> bool
              {
                  return a.second < b.second;
              });


    ss.str("");
    ss.clear();
    for (auto &itr : _modules)
    {
        ss << itr.first->ID << "->";
    }
    s = ss.str();
    SPDLOG_DEBUG("_modules order after sort: {}", s.substr(0, s.length() - 2));

    // Parallel compilation ordering
//    std::vector<int> time(size, 0);
//    for (auto i = make_order.begin(); i != make_order.end(); ++i)
//    {
//        // Walk through the in_edges an calculate the maximum time.
//        if (boost::in_degree(*i, g) > 0)
//        {
//            Graph::in_edge_iterator j, j_end;
//            int maxdist = 0;
//            // Through the order from topological sort, we are sure that every
//            // time we are using here is already initialized.
//            for (boost::tie(j, j_end) = boost::in_edges(*i, g); j != j_end; ++j)
//                maxdist = (std::max)(time[boost::source(*j, g)], maxdist);
//            time[*i] = maxdist + 1;
//        }
//    }
//    boost::graph_traits<Graph>::vertex_iterator i, iend;
//    for (boost::tie(i, iend) = boost::vertices(g); i != iend; ++i)
//    {
//        LOG_DEBUG << "time_slot[" << time[*i] << "] = " << _modules.at(*i).first->ID << std::endl;
//    }



}

void core::_schedule_modules()
{
    //organize modules into sorted parallel data/domain chunks
    size_t chunks = 1; //will be 1 behind actual number as we are using this for an index
    size_t chunk_itr = 0;
    for (auto &itr : _modules)
    {
        SPDLOG_DEBUG( "Chunking module: {}", itr.first->ID);
        //first case, empty list
        if (_chunked_modules.size() == 0)
        {
            _chunked_modules.resize(chunks);
            _chunked_modules.at(0).push_back(itr.first);
        } else
        {
            if (_chunked_modules.at(chunk_itr).at(0)->parallel_type() == itr.first->parallel_type())
            {
                _chunked_modules.at(chunk_itr).push_back(itr.first);
            } else
            {
                chunk_itr++;
                chunks++;
                _chunked_modules.resize(chunks);
                _chunked_modules.at(chunk_itr).push_back(itr.first);
            }
        }
    }

    chunks = 0;
    for (auto &itr : _chunked_modules)
    {
        SPDLOG_DEBUG("Chunk {} {}: ", (itr.at(0)->parallel_type() == module_base::parallel::data ? "data" : "domain"),  chunks);
        for (auto &jtr : itr)
        {
            SPDLOG_DEBUG(jtr->ID);
        }
        chunks++;
    }
}



void core::run()
{
    timer c;

    //setup a XML writer for the PVD paraview format
    pt::ptree pvd;
    pvd.add("VTKFile.<xmlattr>.type", "Collection");
    pvd.add("VTKFile.<xmlattr>.version", "0.1");


    SPDLOG_DEBUG("Loading first timestep's met data");
    // Populate the stations with the first timestep's data.
    // We can do this _once_ without incrementing the internal iterators
    _metdata->next();

    // even though we have already done this once in core::init, if the calling next() above resulted in a load
    // of a new nc on the first time step, e.g., we had to load 2 netcdf to get to the first timestep
    // we need to re build the face-station mapping as  loading a new netcdf will invalidate all our stations
    // TODO: this should probably be moved into metdata, and modules need to know about this to rebuild their
    // interp structs
    if(_metdata->is_multipart_nc() && _metdata->nc_just_loaded())
        populate_face_station_lists();

    SPDLOG_DEBUG("Starting model run");

    c.tic();

    double meantime = 0;
    size_t current_ts = 0;
      size_t max_ts = _metdata->n_timestep();
    bool done = false;

    while (!done)
    {
        boost::posix_time::ptime t;

        _global->timestep_counter++;

        _global->_current_date = _metdata->current_time();

        SPDLOG_DEBUG("Timestep: {}\tstep#{}", boost::posix_time::to_simple_string(_global->posix_time()), current_ts);

        std::stringstream ss;
        ss << _global->posix_time();

        c.tic();
        size_t chunks = 0;
        try
        {
            for (auto &itr : _chunked_modules)
            {

                if (itr.at(0)->parallel_type() == module_base::parallel::data)
                {
#ifdef OMP_SAFE_EXCEPTION
                    ompException e;
#endif
#pragma omp parallel for
                    for (size_t i = 0; i < _mesh->size_local_faces(); i++)
                    {
                        auto face = _mesh->face(i);
                        if (point_mode.enable && face->_debug_name != _outputs[0].name)
                            continue;

                        //module calls
                        for (auto &jtr : itr)
                        {
#ifdef OMP_SAFE_EXCEPTION
                            e.Run(
                                [&]
                                {
#endif
                                    jtr->run(face);
#ifdef OMP_SAFE_EXCEPTION
                                });
#endif
                        }
                    }
#ifdef OMP_SAFE_EXCEPTION
                    e.Rethrow();
#endif

                } else
                {
                    //module calls for domain parallel
                    for (auto &jtr : itr)
                    {
                        jtr->run(_mesh);
                    }
                }

                chunks++;

            }
        }
        catch (exception_base &e)
        {
            SPDLOG_ERROR("Exception at timestep: {}", boost::posix_time::to_simple_string(_global->posix_time()));
            //if we die in a module, try to dump our time series out so we can figure out wtf went wrong
            SPDLOG_ERROR("Exception has occured. Timeseries and meshes WILL BE INCOMPLETE!");
            *_end_ts = _global->posix_time();
            done = true;
            SPDLOG_ERROR(boost::diagnostic_information(e));

        }
        catch(std::exception& e)
        {
            SPDLOG_ERROR("Exception at timestep: {}", boost::posix_time::to_simple_string(_global->posix_time()));
            SPDLOG_ERROR(e.what());
            *_end_ts = _global->posix_time();
            done = true;
            SPDLOG_ERROR(e.what());
        }


        // ensure all ranks agree on this
        bool my_shouldckp = _checkpoint_opts.should_checkpoint(current_ts,
                                               (max_ts-1) == current_ts, // -1 because current_ts is 0 indexed
                                               _hpc_scheduler_info,
                                               _comm_world,
                                               _global->_current_date
                                               );
        bool global_shouldckp=false;
        boost::mpi::all_reduce(_comm_world, my_shouldckp, global_shouldckp, boost::mpi::maximum<bool>());
        // save the current state
        if(global_shouldckp) //
        {
            SPDLOG_DEBUG("Checkpointing...");

            netcdf savestate; //file to save to when checkpointing.

            auto timestamp = _global->posix_time() + boost::posix_time::seconds(_global->_dt);
            //also write it out in seconds because netcdf is struggling with the string
            unsigned long long int ts_sec = _global->posix_time_int()+_global->_dt;

            auto timestr = boost::posix_time::to_iso_string(timestamp); // start from current TS + dt

            auto dirpath = _checkpoint_opts.ckpt_path / timestr;
            if (_comm_world.rank() == 0)
            {
                boost::filesystem::create_directories(dirpath);
            }

            _comm_world.barrier();

            //this parses both the input and the output paths for the checkpoint.
            auto fname = ("chkp"+timestr + "_" + std::to_string(_comm_world.rank()) + ".nc");
            auto f = dirpath / fname;
            savestate.create( f.string());

            c.tic();
            for (auto &itr : _chunked_modules)
            {
                //module calls
                for (auto &jtr : itr)
                {
                    jtr->checkpoint(_mesh, savestate);
                }
            }

            auto& ids = _mesh->get_global_IDs();
            savestate.create_variable1D("global_id",ids.size());

            for (size_t i = 0; i < ids.size(); i++)
            {
                savestate.put_var1D("global_id", i, ids[i]);
            }

            savestate.get_ncfile().putAtt("restart_time",boost::posix_time::to_simple_string(timestamp));
            savestate.get_ncfile().putAtt("restart_time_sec", netCDF::ncUint64, ts_sec);
            savestate.get_ncfile().putAtt("timestep_counter", netCDF::ncUint64, static_cast<unsigned long long>(_global->timestep_counter));

            pt::ptree tree;

            int nranks = _comm_world.size();


            tree.put("ranks", nranks);
            tree.put("restart_time_sec", ts_sec);
            tree.put("startdate", timestr);
            tree.put("timestep_counter", _global->timestep_counter);

            pt::ptree files;

            pt::ptree tmp_files;
            for (int i = 0; i < nranks; ++i)
            {
                pt::ptree s;

                s.put("", timestr +"/" + "chkp"+timestr + "_" + std::to_string(i) + ".nc");
                tmp_files.push_back(std::make_pair("", s));
            }
            tree.add_child("files", tmp_files);


            if(_comm_world.rank() == 0)
            {
                pt::write_json(
                    (_checkpoint_opts.ckpt_path / ("checkpoint_" + timestr + ".np" + std::to_string(nranks) + ".json")).string(),
                    tree);
            }

            SPDLOG_DEBUG("Done checkpoint [ {} s]", c.toc<s>());

            // if we checkpointed because we are out of time, we need to stop the simulation
            if(_checkpoint_opts.checkpoint_request_terminate)
            {
                // we bailed early because of wall clock, so this is not a clean exit
                // done = true; clean_exit = false;

                // ensure everyone agrees we should be done = true
                boost::mpi::all_reduce(_comm_world, true, done, boost::mpi::maximum<bool>());
                // ensure everyone aggress we should be clean_exit = false
                boost::mpi::all_reduce(_comm_world, false, clean_exit, boost::mpi::minimum<bool>());
            }
        }

        for (auto &itr : _outputs)
        {
            if (itr.type == output_info::output_type::mesh)
            {
                // check if we should output or not
                bool do_output = itr.should_output(max_ts, current_ts, _global->_current_date);

                if(do_output)
                {
                    _mesh->set_output_parameters(itr.output_parameters);

                    std::vector<std::string> output;
                    output.assign(itr.variables.begin(),itr.variables.end()); //convert to list to match internal lists

                    if (itr.mesh_output_formats == output_info::mesh_outputs::vtu)
                    {
                        std::string base_name = itr.fname + std::to_string(_global->posix_time_int());
                        boost::filesystem::path p(base_name);
                        _mesh->update_vtk_data(output); //update the internal vtk mesh

                        // this really only works if we let rank0 handle the io.
                        // If we let each process do it, they walk all over each other's output

                        if(_comm_world.rank() == 0)
                        {
                            for(int rank = 0; rank < _comm_world.size(); rank++)
                            {

                                // write paths that are relative to the pvd file
                                boost::filesystem::path vtu_path(output_folder_path.string() + "/vtu/" + p.filename().string()+"_"+std::to_string(rank) + ".vtu");
                                pt::ptree &dataset = pvd.add("VTKFile.Collection.DataSet", "");
                                dataset.add("<xmlattr>.timestep", _global->posix_time_int());
                                dataset.add("<xmlattr>.group", "");
                                dataset.add("<xmlattr>.part", rank);
                                dataset.add("<xmlattr>.file", boost::filesystem::relative(vtu_path, output_folder_path).string());

                            }
                        }

                        //because a full path can be provided for the base_name, we need to strip this off
                        //to make it a relative path in the xml file.
                        _mesh->write_vtu(base_name + "_"+std::to_string(_comm_world.rank() )+ ".vtu");

                    }
                    else if (itr.mesh_output_formats == output_info::mesh_outputs::ugrid)
                    {
                        // first, check if we need a new ugrid file
                        bool new_ugrid = itr.should_rotate(max_ts, current_ts, _global->_current_date);

                        auto& writer = boost::get<boost::shared_ptr<ugrid_writer>>(itr.writer);
                        if (new_ugrid)
                        {
                            writer->close_ugrid();
                        }

                        writer->write_ugrid({itr.variables.begin(), itr.variables.end()} );
                    }
                }
            }
        }

        //If we are output a timeseries at specific triangles, we do that here
        //Each output knows what face it corresponds to
        for (auto &itr : _outputs)
        {
            //only update the full timeseries
            if (itr.type == output_info::output_type::time_series)
            {
                for (const auto& v : _provided_var_module)
                {
                    auto data = (*itr.face)[v];
                    itr.ts.at(v, current_ts) = data;
                }
            }
        }

        if(!_metdata->next())
        {
            done = true;
        }
        else
        {
            // loading a new netcdf will invalidate all our stations, so we need to rebuild the list of stations
            if(_metdata->is_multipart_nc() && _metdata->nc_just_loaded())
                populate_face_station_lists();
        }


        auto timestep = c.toc<ms>();
        meantime += timestep;

        current_ts++;


        double mt = meantime / current_ts;
        bool ms = true;
        if (mt > 1000)
        {
            mt /= 1000.;
            ms = false;
        }

        std::string s = std::to_string(std::lround(mt)) + (ms == true ? " ms" : "s");


        //we need it in seconds now
        if (ms)
        {
            mt /= 1000.0;
        }

        boost::posix_time::ptime pt(boost::posix_time::second_clock::local_time());
        pt = pt + boost::posix_time::seconds(size_t(mt) * (max_ts - current_ts));

        SPDLOG_DEBUG("Took {}s. Avg duration {} \tEstimated completion: {}", std::lround(timestep/1000), s, boost::posix_time::to_simple_string(pt));
        _global->first_time_step = false;


    }
    double elapsed = c.toc<s>();
    SPDLOG_DEBUG("Total runtime was {}s", elapsed);



    std::string base_name="";

    for (auto &itr : _outputs)
    {
        if (itr.mesh_output_formats == output_info::mesh_outputs::vtu)
        {

#ifdef USE_MPI
            if(_comm_world.rank() == 0)
            {
#endif

                // output the pvd one level higher in the main outdir than we have previously
                boost::filesystem::path path(itr.fname + ".pvd");
                pt::write_xml( (output_folder_path.string() / path.filename()).string(),
                              pvd, std::locale(), pt::xml_writer_settings<std::string>(' ', 4));
                break;

#ifdef USE_MPI
            }
#endif
        }
    }

    for (auto &itr : _outputs)
    {
        //save the full timeseries
        if (itr.type == output_info::output_type::time_series)
        {
            itr.ts.subset(*_start_ts,*_end_ts); // in the event of an exception, _end_ts will be reset to have the esception timestep so-as to no write massive amounts of nan values
            itr.ts.to_file(itr.fname);
        }

        if (itr.mesh_output_formats == output_info::mesh_outputs::ugrid)
        {
            boost::get<boost::shared_ptr<ugrid_writer>>(itr.writer)->close_ugrid();
        }
    }

    if(_notification_script != "")
    {
        SPDLOG_DEBUG("Calling notification script");
        int ierr = std::system(_notification_script.c_str()); CHK_SYSTEM_ERR(ierr);
    }
}

void core::end(const bool abort)
{

    // Write the sentinel file IFF there is a clean exit
    // recall that if the checkpointing system detects an about-to-expire wallclock and terminates early,
    // this doesn't count as a clean exit.
    if(!abort && clean_exit)
    {
        if(_comm_world.rank() == 0)
        {
            std::ofstream((output_folder_path / "clean_exit").string()).close();
        }
    }
    SPDLOG_DEBUG("Finished cleaning up");


    if(abort)
    {

        SPDLOG_ERROR("An exception has occurred, requesting MPI Abort! Log files of failed will be in <output>/error_logs");

        // best-effort copy of the log before aborting nukes the process
        // mpi abort will trigger a shutdown and avoid the descructors so we have to process the log movement here
        try
        {
            // flush anything pending before we move the log files
            try
            {
                spdlog::default_logger()->flush();
            }
            catch (...) {}

            std::string job_name = _hpc_scheduler_info.job_name;
            if (job_name.empty())
            {
                job_name = "unknown.jobid";
            }
            const auto err_path = output_folder_path / ("error_logs_" + job_name);
            const auto log_path = output_folder_path / ("logs_" + job_name);

            try { boost::filesystem::create_directories(err_path); } catch (...) {}
            try { boost::filesystem::create_directories(log_path); } catch (...) {}

            try
            {
                boost::filesystem::copy_file(
                    log_file_path,
                    err_path / log_file_path.filename(),
                    boost::filesystem::copy_options::overwrite_existing);
            }
            catch (...) {}

            try
            {
                boost::filesystem::rename(
                    log_file_path,
                    log_path / log_file_path.filename());
            }
            catch (...) {}
        }
        catch (...)
        {
            // swallow all errors; we are aborting regardless
        }

        _mpi_env.abort(-1);
    }

}

bool core::check_is_geographic(const std::string& path)
{
    auto mesh_file_extension = boost::filesystem::path(path).extension();

    std::string mesh_path = path;
    bool is_geographic = false;

    // we have a partition so load just the first part of the partition
    if(mesh_file_extension == ".partition")
    {
        try
        {
            auto json = read_json(path);
            for(auto &itr : json.get_child("meshes"))
            {
                auto mesh_file_paths = boost::filesystem::path(itr.second.data());
                if(mesh_file_paths.is_relative())
                    mesh_path = (cwd_dir / mesh_file_paths).string();
                else
                    mesh_path = mesh_file_paths.string();
                mesh_file_extension = boost::filesystem::path(mesh_path).extension();
                break;
            }
        }
        catch(pt::ptree_bad_path &e)
        {
            CHM_THROW_EXCEPTION(config_error, "A mesh is required and was not found");
        }
    }

    if(mesh_file_extension == ".h5")
    {
        hsize_t geographic_dims = 1;
        H5::DataSpace dataspace(1, &geographic_dims);

        H5File  file(mesh_path, H5F_ACC_RDONLY);
        H5::Attribute attribute = file.openAttribute("/mesh/is_geographic");
        attribute.read(PredType::NATIVE_HBOOL, &is_geographic);
        file.close();
    }
    else
    {
        auto mesh = read_json(path);
        if( mesh.get<int>("mesh.is_geographic") == 1)
        {
            is_geographic = true;
        }
        else
        {
            is_geographic = false;
        }
    }

    return is_geographic;
}

void core::populate_face_station_lists()
{

    SPDLOG_DEBUG("Populating each face's station list");

    for (size_t i = 0; i < _mesh->size_local_faces(); i++)
    {
        auto f = _mesh->face(i);
        f->stations().clear();

        auto stations = _metdata->get_stations(f->get_x(), f->get_y());
        f->stations().insert(std::end(f->stations()), std::begin(stations), std::end(stations));

        f->nearest_station() = _metdata->nearest_station(f->get_x(), f->get_y()).at(0);
    }

}

std::vector< std::pair<module,size_t> >& core::get_active_module_list()
{
    return _modules;
}
