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
    BOOST_LOG_FUNCTION();

    _start_ts = nullptr;
    _end_ts = nullptr;
    _interpolation_method = interp_alg::tpspline;

    //default logging level
    _log_level = debug;
    _output_station_ptv = true;
    _use_netcdf=false;
    _load_from_checkpoint=false;
    _do_checkpoint=false;
    _metdata= nullptr;
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

}

void core::config_options( pt::ptree &value)
{
    LOG_DEBUG << "Found options section";

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

    LOG_DEBUG << "Setting log severity to " << _log_level;

    _log_sink->set_filter(
            severity >= _log_level
    );

    _cout_log_sink->set_filter(
            severity >= _log_level
    );


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
        LOG_WARNING << "Unknown interpolant selected, defaulting to spline";
    }

    // custom start time
    boost::optional<std::string> start = value.get_optional<std::string>("startdate");
    if (start)
    {
        _start_ts = new boost::posix_time::ptime(boost::posix_time::from_iso_string(*start));
        LOG_DEBUG << "User-specified startdate: " << *_start_ts;
    }

    // custom start time
    boost::optional<std::string> end = value.get_optional<std::string>("enddate");
    if (end)
    {
        _end_ts = new boost::posix_time::ptime(boost::posix_time::from_iso_string(*end));
        LOG_DEBUG << "User-specified endate: " << *_end_ts;
    }


    // point mode options
    try
    {
        auto pm = value.get_child("point_mode");

        point_mode.enable = true;
        point_mode.forcing = pm.get<std::string>("forcing");
        point_mode.output  = pm.get<std::string>("output");
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

    auto radius = value.get_optional<double>("station_search_radius");
    auto N = value.get_optional<double>("station_N_nearest");

    if(radius && N)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info("Cannot have both station_search_radius and station_N_nearest set."));
    }

    if(radius)
    {
        _metdata->get_stations = boost::bind( &metdata::get_stations_in_radius,_metdata,boost::placeholders::_1,boost::placeholders::_2, *radius);
    }
    else
    {
        int n = 0 ; // Number of stations to interp
        if(N) // If user specified N in config
        {
            n = *N;
        }
        else
        { // N not specified used defaults
            if (ia == "nearest")
            {
                n = 1;
                LOG_WARNING << "Using N=1 nearest stations as default.";
            } else
            {
                n = 5;
                LOG_WARNING << "Using N=5 nearest stations as default.";
            }
        }

        if( (n < 2) && (ia != "nearest")) // Required more than 1 station if using spline or idw
        {
            BOOST_THROW_EXCEPTION(config_error() << errstr_info("station_N_nearest must be >= 2 if spline or idw is used. N = " + std::to_string(n)));
        }

        _metdata->get_stations = boost::bind( &metdata::nearest_station,_metdata,boost::placeholders::_1,boost::placeholders::_2, n);
    }


}

void core::config_module_overrides( pt::ptree &value)
{
    _find_and_insert_subjson(value);
    LOG_DEBUG << "Found dependency override section";
    for (auto &itr : value)
    {
        std::string A = itr.first.data();
        std::string B = itr.second.data();

        LOG_WARNING << "Removing depency of " << B << " from " << A;
        _overrides.push_back(std::make_pair(A, B));

    }
}

void core::config_modules(pt::ptree &value, const pt::ptree &config, std::vector<std::string> remove,
                          std::vector<std::string> add)
{
    LOG_DEBUG << "Found modules section";
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
            LOG_DEBUG << "Removed module " << module;
        }

    }

    //straight up add anything
    for (auto &itr : add)
    {
        modules.insert(itr);
        LOG_DEBUG << "Inserted module " << itr << " from cmdl";
    }


    for (auto &itr : modules)
    {

        std::string module_name = itr;
        LOG_DEBUG << "Module ID=" << module_name;


        //try grabbing a config for this module, empty string default
        pt::ptree cfg;
        try
        {
            cfg = config.get_child(module_name);
        } catch (pt::ptree_bad_path &e)
        {
            LOG_DEBUG << "No config for " << module_name;
        }

        boost::shared_ptr<module_base> module = module_factory::create(module_name,cfg);
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
        BOOST_THROW_EXCEPTION(no_modules_defined() << errstr_info("No modules defined. Aborting"));
    }
}

void core::config_checkpoint( pt::ptree& value)
{
    LOG_DEBUG << "Found checkpoint section";
    _find_and_insert_subjson(value);

    _do_checkpoint = value.get("save_checkpoint",false);

    //this uses the output path defined in config_output, so this must run after that.
    if(_do_checkpoint)
    {
        auto dir = "checkpoint";

        boost::filesystem::path ckpt_path;

        // Create empty folder
        ckpt_path = o_path / dir;
        boost::filesystem::create_directories(ckpt_path);

        //this parses both the input and the output paths for the checkpoint.
        auto f = ckpt_path / "checkpoint.nc";
        _savestate.create( f.string());

        _checkpoint_feq = value.get("frequency",1);
        LOG_DEBUG << "Checkpointing every " << _checkpoint_feq << " timesteps.";

    }

    auto file = value.get_optional<std::string>("load_checkpoint_path");

    if (file)
    {
        _load_from_checkpoint = true;

        boost::filesystem::path ckpt_path;
        ckpt_path = cwd_dir / *file;

        _checkpoint_file = ckpt_path.string();

        _in_savestate.open(_checkpoint_file);
    }



}
void core::config_forcing(pt::ptree &value)
{
    LOG_DEBUG << "Found forcing section";
    LOG_DEBUG << "Reading meta data from config";

    //positive offset going west. So the normal UTC-6 would be UTC_offset:6
    _global->_utc_offset = value.get("UTC_offset",0);
    LOG_DEBUG << "Applying UTC offset to ALL forcing files. UTC+" << std::to_string(_global->_utc_offset);

    _find_and_insert_subjson(value);

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

        // this delegates all filter responsibility to metdata from now on
        _metdata->load_from_netcdf(file, netcdf_filters);
        nstations = _metdata->nstations();
    } else
    {
        std::vector<metdata::ascii_metdata> ascii_data;

        for (auto &itr : value)
        {
            if(itr.first != "UTC_offset")
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

    LOG_DEBUG << "Found # stations = " <<  nstations;
    if(nstations == 0)
    {
        CHM_THROW_EXCEPTION(forcing_error,"No input forcing files found!");
    }


    auto f = o_path / "stations.vtp";
    _metdata->write_stations_to_ptv(f.string());

    LOG_DEBUG << "Finished reading stations. Took " << c.toc<s>() << "s";

}
void core::determine_startend_ts_forcing()
{

    if (_metdata->nstations() == 0)
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("no stations"));

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
    if(_load_from_checkpoint)
    {
        size_t t = 0;
        _in_savestate.get_ncfile().getAtt("restart_time_sec").getValues(&t);
        _start_ts = new boost::posix_time::ptime(boost::posix_time::from_time_t(t));

        LOG_WARNING << "Loading from checkpoint. Overriding start time to match. New start time = " << *_start_ts;
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
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info(ss.str()));
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
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info(ss.str()));
    }


    if( *_start_ts > *_end_ts)
    {
        std::stringstream ss;
        ss << "Start time is after endtime. ";
        ss << " Start date: " << *_start_ts;
        ss << " End date: " << *_end_ts;
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info(ss.str()));
    }


    //figure out what our timestepping is. Needs to happen before we subset as we may end up with only 1 timestep
    _global->_dt = _metdata->dt_seconds();
    LOG_DEBUG << "model dt = " << _global->dt() << " (s)";


    LOG_DEBUG << "Subsetting station timeseries";
    _metdata->subset(*_start_ts, *_end_ts);


    //ensure all the stations have the same start and end times
    // per-timestep agreement happens during runtime.
    _metdata->check_ts_consistency();
}
void core::config_parameters(pt::ptree &value)
{
    LOG_DEBUG << "Found parameter mapping section";

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

void core::config_meshes( pt::ptree &value)
{
    LOG_DEBUG << "Found meshes sections.";
#ifdef MATLAB
    _mesh = boost::make_shared<triangulation>(_engine);
#else
    _mesh = boost::make_shared<triangulation>();
#endif

    _mesh->_global = _global;

    _find_and_insert_subjson(value);

    std::string mesh_path = value.get<std::string>("mesh");
    LOG_DEBUG << "Found mesh:" << mesh_path;
    mesh_path = (cwd_dir / mesh_path).string();

    auto mesh_file_extension = boost::filesystem::path(mesh_path).extension().string();

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
                LOG_DEBUG << "Found parameter file: " << param_file;
                param_file_paths.push_back( (cwd_dir / itr.second.data()).string() );
            }
        }
        catch(pt::ptree_bad_path &e)
        {
            LOG_DEBUG << "No addtional parameters found in mesh section.";
        }

        // Paths for initial condition files
        try
        {
            for(auto &itr : value.get_child("initial_conditions"))
            {
                auto ic_file = itr.second.data();
                LOG_DEBUG << "Found initial condition file: " << ic_file;
                initial_condition_file_paths.push_back( (cwd_dir / itr.second.data()).string() );
            }
        }
        catch(pt::ptree_bad_path &e)
        {
            LOG_DEBUG << "No addtional initial conditions found in mesh section.";
        }
    }


    // Ensure all files are HDF5 in multiprocess MPI runs
#ifdef USE_MPI
    if ( _comm_world.size() > 1 )
    {
        bool h5_or_part = mesh_file_extension != ".h5" || mesh_file_extension != ".partition";
        if(!h5_or_part)
        {
            BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("MPI multiprocess run requires hdf5 mesh.\n\n    Run the serial hdf5 conversion tool\n\n"));
        }

        // only check the params and ics if we aren't using a parition file
        if(mesh_file_extension != ".partition")
        {
            for(const auto& it : param_file_paths)
            {
                auto extension = boost::filesystem::path(it).extension();
                if(extension != ".h5")
                {
                    BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("MPI multiprocess run requires hdf5 parameter files: " + it + "\n\n    Run the serial hdf5 conversion tool\n\n"));
                }
            }
            for(const auto& it : initial_condition_file_paths)
            {
                auto extension = boost::filesystem::path(it).extension();
                if(extension != ".h5")
                {
                    BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("MPI multiprocess run requires hdf5 initial condition files: " + it + "\n\n    Run the serial hdf5 conversion tool\n\n"));
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
    bool is_geographic = check_is_geographic(mesh_path);

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
        _mesh->from_hdf5(mesh_path, param_file_paths, initial_condition_file_paths);
    }
    else if(mesh_file_extension == ".partition")
    {
#ifndef USE_MPI
        CHM_THROW_EXCEPTION(config_error, "Using a partitioned mesh requires enabling MPI during CHM configure and build.");
#endif

        _mesh->from_partitioned_hdf5(mesh_path);



    }
    else  // Assume anything that is NOT h5 is json
    {
        auto mesh = read_json(mesh_path);

        //mesh will have been loaded by the geographic check so don't re load it here
      bool triarea_found = false;

      // Parameter files
      for (auto param_file : param_file_paths)
      {
	    pt::ptree param_json = read_json(param_file);

            for(auto& ktr : param_json)
            {
                //use put to ensure there are no duplciate parameters...
                std::string key = ktr.first.data();
                mesh.put_child( "parameters." + key ,ktr.second);

                if( key == "area")
                    triarea_found = true;

                LOG_DEBUG << "Inserted parameter " << ktr.first.data() << " into the config tree.";
	    }
      }

        if(is_geographic && !triarea_found)
        {
            BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Geographic meshes require the triangle area be present in a .param file. Please include this."));
        }

      // Initial condition files
      for(auto ic_file : initial_condition_file_paths)
      {
            pt::ptree ic_json = read_json(ic_file);

            for(auto& ktr : ic_json)
            {
                //use put to ensure there are no duplciate parameters...
                std::string key = ktr.first.data();
                mesh.put_child( "initial_conditions." + key ,ktr.second);
                LOG_DEBUG << "Inserted initial condition " << ktr.first.data() << " into the config tree.";
            }
      }

      _mesh->from_json(mesh);

    }

    if (_mesh->size_faces() == 0)
        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Mesh size = 0!"));

}

void core::config_matlab( pt::ptree &value)
{
    LOG_DEBUG << "Found matlab section";
#ifdef MATLAB
    //loop over the list of matlab options
    for (auto& jtr : value.get_obj())
    {
        const json_spirit::Pair& pair = jtr;
        const std::string& name = pair.name_;
        const json_spirit::Value& value = pair.value_;

        if (name == "mfile_paths")
        {
            LOG_DEBUG << "Found " << name;
            for (auto& ktr : value.get_obj()) //loop over all the paths
            {
                const json_spirit::Pair& pair = ktr;
                //                const std::string& name = pair.name_;
                const json_spirit::Value& value = pair.value_;


                _engine->add_dir_to_path(value.get_str());

            }
        }
    }
#endif
}

void core::config_output(pt::ptree &value)
{
    LOG_DEBUG << "Found output section";
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkStringArray> labels = vtkSmartPointer<vtkStringArray>::New();
    labels->SetName("Point output name");
    _output_station_ptv = true;

    size_t ID = 0;
    auto pts_dir = "points";
    auto msh_dir = "meshes";
    boost::filesystem::path pts_path;
    boost::filesystem::path msh_path;


    auto output_dir = value.get<std::string>("output_dir","output");
    o_path = cwd_dir / output_dir;
    boost::filesystem::create_directories(o_path);

    // Create empty folders /points/ and /meshes/
    pts_path = o_path / pts_dir;
    msh_path = o_path / msh_dir;
    boost::filesystem::create_directories(pts_path);
    boost::filesystem::create_directories(msh_path);

    _find_and_insert_subjson(value);

    for (auto &itr : value)
    {
        output_info out;
        std::string out_type = itr.first.data();
        if (out_type == "output_dir") //skip this
        {
            continue;
        }

        if ((out_type != "mesh"))  // anything else *should* be a time series*......
        {
            out.type = output_info::time_series;
            out.name = out_type;

            std::string fname = "";
            try
            {
                fname = itr.second.get<std::string>("file");
            }
            catch(const pt::ptree_error &e)
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Missing output filename for " + out.name));
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
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Output point " + out.name + " is missing latitude and/or longitude."));
            }


            if( (out.latitude > 90 || out.latitude < -90) ||
                (out.longitude > 180 || out.longitude < -180) )
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Output " + out.name + " coordinate is invalid."));
            }

            //project mesh, need to convert the input lat/long into the coordinate system our mesh is in
            if(!_mesh->is_geographic())
            {

                OGRSpatialReference insrs;
                insrs.SetWellKnownGeogCS("WGS84");

                OGRSpatialReference outsrs;
                outsrs.importFromProj4(_mesh->proj4().c_str());


                OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&insrs, &outsrs);

                //CRS created with the “EPSG:4326” or “WGS84” strings use the latitude first, longitude second axis order.
                if(!coordTrans->Transform(1, &out.latitude, &out.longitude))
                {
                    BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Output=" + out.name + ": unable to convert coordinates to mesh format."));
                }

                OGRCoordinateTransformation::DestroyCT(coordTrans);
            }

            out.face = _mesh->locate_face(out.longitude, out.latitude);

            if(out.face != nullptr)
            {
                // In MPI mode, we might not thave the output triangle on this node, so we just have to fail gracefully and hope that the other nodes have this.
                // in the future we need a comms here to check if the nodes successfully figured out the output points
                // for now, if we don't have it, print an error but keep going

                //set the point to be the center of the triangle that the output point lies on
                points->InsertNextPoint(out.face->get_x(), out.face->get_y(), out.face->get_z());
                labels->InsertNextValue(out.name);

                LOG_DEBUG << "Triangle geometry for output triangle = " << out_type << " slope: "
                          << out.face->slope() * 180. / M_PI << " aspect:" << out.face->aspect() * 180. / M_PI;

                out.face->_debug_name = out.name; //out_type holds the station name
                out.face->_debug_ID = ID;
                ++ID;
            }
        }
        else if (out_type == "mesh")
        {
            out.type = output_info::mesh;

            auto fname = itr.second.get<std::string>("base_name","output");
            auto f = msh_path / fname;
            boost::filesystem::create_directories(f.parent_path());
            out.fname = f.string();

            _mesh->write_param_to_vtu( itr.second.get("write_parameters",true) ) ;

	    // Set option for writing ghost neighbor data, defaults to not
            _mesh->write_ghost_neighbors_to_vtu( itr.second.get("write_ghost_neighbors",false) ) ;

            try
            {
                for (auto &jtr: itr.second.get_child("variables"))
                {
                    out.variables.insert(jtr.second.data());
                }
            }
            catch (pt::ptree_bad_path &e)
            {
                LOG_DEBUG << "Writing all variables to output mesh";
            }

            out.frequency = itr.second.get("frequency",1); //defaults to every timestep
            LOG_DEBUG << "Output every " << out.frequency <<" timesteps.";

            out.only_last_n = itr.second.get("only_last_n",-1); //defaults to keep everything
            LOG_DEBUG << "Output only on last n =  " << out.only_last_n <<" timesteps.";

            if(out.frequency > 1 && out.only_last_n != -1)
            {
                LOG_WARNING << "Only only_last_n output option will be used";
            }

            out.mesh_output_formats.push_back(output_info::mesh_outputs::vtu);

        } else
        {
            LOG_WARNING << "Unknown output type: " << itr.second.data();
        }

        //ensure we aren't adding timeseries outputs with invalid faces bc of MPI
        if(  out.type == output_info::time_series &&
             out.face == nullptr)
        {
            // we will pass for now, but ultimiately we should have done an MPI comms and
            // check if we are missing an output
#ifndef USE_MPI
            BOOST_THROW_EXCEPTION(config_error() <<
                                                     errstr_info(
                                                             "Requested an output point that is not in the triangulation domain. Pt:"
                                                             + std::to_string(out.longitude) + "," +
                                                             std::to_string(out.latitude) + " name: " + out.name));
#else
            LOG_WARNING << "In MPI mode there is currently no check if all the nodes correctly find the output triangle. "
                           "If you are missing output, ensure that all the output points are within the domain.";
#endif

        }else
        {
            _outputs.push_back(out);
        }
    }
#ifdef USE_MPI
    LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has #ouput points = " << _outputs.size();
#endif
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(labels);
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    //output this to the same folder as the points are written out to
    std::string rank = "";
#ifdef USE_MPI
    rank = "."+ std::to_string(_comm_world.rank());
#endif
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
    BOOST_THROW_EXCEPTION(no_modules_defined() << errstr_info("Global section is now deprecated. This does nothing, please remove from config file."));
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
                BOOST_THROW_EXCEPTION(io_error() << errstr_info("Config value of " + itr + " is invalid."));

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
        LOG_ERROR << "Configuration file required.";
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
    BOOST_LOG_FUNCTION();


    //default logging level
    _log_level = debug;


    _cout_log_sink = boost::make_shared<text_sink>();

    text_sink::locked_backend_ptr pBackend_cout = _cout_log_sink->locked_backend();

#if (BOOST_VERSION / 100 % 1000) < 56
    boost::shared_ptr< std::ostream > pStream(&std::clog,  logging::empty_deleter());
#else
    boost::shared_ptr<std::ostream> pStream(&std::cout, boost::null_deleter()); //clog
#endif
    pBackend_cout->add_stream(pStream);

    _cout_log_sink->set_formatter
            (
                    expr::format("[%1%]: %2%")
                    % expr::attr<log_level>("Severity")
                    % expr::smessage
            );
    _cout_log_sink->set_filter(
            severity >= _log_level
    );
    logging::core::get()->add_sink(_cout_log_sink);



    auto log_start_time = boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());

    //get any command line options
    auto cmdl_options = config_cmdl_options(argc, argv);

    //see if we are getting passed a path for our config file. If so, we need to append this path to
    //any config files we are about to read in the main config file.

    //do this here though as we need cwd_dir for the log output path
    boost::filesystem::path path(cmdl_options.get<0>());
    cwd_dir = boost::filesystem::current_path();

    LOG_DEBUG << "Current working directory: " << cwd_dir.string();


    std::string log_dir = "log";
    auto log_path = cwd_dir / log_dir;

    //output a unique logfile for each mpi rank
    std::string rank = "";
#ifdef USE_MPI
   rank = "."+std::to_string(_comm_world.rank());
#endif


    std::string log_name = "CHM_" + log_start_time + rank + ".log";

    boost::filesystem::create_directories(log_path);

    //fully qualified path + fname to use in the ofstream calls
    if(cmdl_options.get<5 >()) //-legacy-log true option
    {

        log_name =  "CHM.log";
        log_file_path = cwd_dir / log_name;
    }
    else {
        log_file_path = log_path / log_name;
        log_name =  log_file_path.string();
    }


    _log_sink = boost::make_shared<text_sink>();
    text_sink::locked_backend_ptr pBackend_file = _log_sink->locked_backend();
    boost::shared_ptr<std::ofstream> pStream2(new std::ofstream(log_name));

    if (!pStream2->is_open())
    {
        BOOST_THROW_EXCEPTION(file_write_error()
                              << boost::errinfo_errno(errno)
                              << boost::errinfo_file_name(log_name)
        );
    }

    pBackend_file->add_stream(pStream2);
    pBackend_file->auto_flush (true);

    _log_sink->set_formatter
            (
                    expr::format("%1% [%2%]: %3%")
//                    % expr::attr<boost::posix_time::ptime>("TimeStamp")
                    % expr::format_named_scope("Scope",
                                               keywords::format = "%n:%l",
                                               keywords::iteration = expr::reverse,
                                               keywords::depth = 2)
                    % expr::attr<log_level>("Severity")
                    % expr::smessage
            );
    _log_sink->set_filter(
            severity >= verbose
    );

//    logging::core::get()->add_global_attribute("TimeStamp", attrs::local_clock());
    logging::core::get()->add_global_attribute("Scope", attrs::named_scope());



    logging::core::get()->add_sink(_log_sink);

    LOG_DEBUG << "Logger initialized. Writing to cout and " + log_name;
    LOG_DEBUG << version;



#ifdef NOMATLAB
    _engine = boost::make_shared<maw::matlab_engine>();


    _engine->start();
    _engine->set_working_dir();

    LOG_DEBUG << "Matlab engine started";
#endif

    _global = boost::make_shared<global>();

    // This needs to be set so that underflows in gsl math
    // computations are not treated as errors. i.e.,
    //    gsl: expint.c:363: ERROR: underflow
    //    Default GSL error handler invoked.
    // is not what we want.
    gsl_set_error_handler_off();

#ifdef USE_MPI
    LOG_DEBUG << "Built with MPI support,    #processes = " << _comm_world.size();
#endif

#ifdef _OPENMP
    LOG_DEBUG << "Built with OpenMP support, #threads   = " << omp_get_max_threads();
#endif

    // Check if we are running under slurm
    const char* SLURM_JOB_ID = std::getenv("SLURM_JOB_ID");
    if(SLURM_JOB_ID)
    {
        const char* SLURM_TASK_PID = std::getenv("SLURM_TASK_PID"); //The process ID of the task being started.
        const char* SLURM_PROCID = std::getenv("SLURM_PROCID"); // The MPI rank (or relative process ID) of the current process

        LOG_DEBUG << "Detected running under SLURM as jobid " << SLURM_JOB_ID;
        LOG_DEBUG << "SLURM_TASK_PID = " << SLURM_TASK_PID;
        LOG_DEBUG << "SLURM_PROCID = " << SLURM_PROCID;
    }



    pt::ptree cfg;
    try
    {
        LOG_DEBUG << "Reading configuration file " << cmdl_options.get<0>();
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
        LOG_DEBUG << "Optional section Module config not found";
    }
    catch (pt::json_parser_error &e)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                "Error reading file: " + cmdl_options.get<0>() + " on line: " + std::to_string(e.line()) + " with error: " +
                e.message()));
    }



    //now apply any override or extra configuration parameters from the command line
    for (auto &itr : cmdl_options.get<1>())
    {
        //check if we are overwriting something
        try
        {
            auto value = cfg.get<std::string>(itr.first);
            LOG_WARNING << "Overwriting " << itr.first << "=" << value << " with " << itr.first << "=" <<
                        itr.second;

        } catch (pt::ptree_bad_path &e)
        {
            LOG_DEBUG << "Inserting new config " << itr.first << "=" << itr.second;
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


            LOG_DEBUG << "Removing " << itr;

        } catch (pt::ptree_bad_path &e)
        {
            LOG_DEBUG << "No value " << itr << " to remove";
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
    config_meshes(cfg.get_child("meshes")); // this must come before forcing, as meshes initializes the required distance functions based on geographic/utm meshes

    // This needs to be initialized with the mesh prior to the forcing and output being dealt with.
    // met data needs to know about the meshes' coordinate system. Probably worth pulling this apart further
    _metdata = std::make_shared<metdata>(_mesh->proj4());

    //output should come before forcing, controls if we should output the vtp file of station locations
    try
    {
        config_output(cfg.get_child("output"));
    } catch (pt::ptree_bad_path &e)
    {
        _output_station_ptv = false;
        LOG_DEBUG << "Optional section Output not found";
    }

    config_forcing(cfg.get_child("forcing"));

    /*
     * We can expect the following sections to be optional.
     */
    try
    {
        config_options(cfg.get_child("option"));
    } catch (pt::ptree_bad_path &e)
    {
        LOG_DEBUG << "Optional section option not found";
    }


    // INSERT STATION TRIMMING HERE (after options for interpolation stuff has occurred)
    populate_face_station_lists();
    populate_distributed_station_lists();


    boost::filesystem::path full_path(boost::filesystem::current_path());

    try
    {
        config_global(cfg.get_child("global"));
    } catch (pt::ptree_bad_path &e)
    {
        LOG_DEBUG << "Optional section Global not found";
    }

    try
    {
        config_module_overrides(cfg.get_child("remove_depency"));
    } catch (pt::ptree_bad_path &e)
    {
        LOG_DEBUG << "Optional section remove_depency not found";
    }

    try
    {
        config_parameters(cfg.get_child("parameter_mapping"));
    } catch (pt::ptree_bad_path &e)
    {
        LOG_DEBUG << "Optional section parameter mapping not found";
    }

    try
    {
        config_checkpoint(cfg.get_child("checkpoint"));
    } catch (pt::ptree_bad_path &e)
    {
        LOG_DEBUG << "Optional section checkpoint not found";
    }

//#ifdef NOMATLAB
//            config_matlab(value);
//#endif

    //if we are to run in point mode, we need to remove all the input and outputs that aren't associated with point mode
    // we do this so the user doesn't have to modify a lot of the config file when going between point and dist mode.

    if(point_mode.enable)
    {
        LOG_INFO << "Running in point mode";

        std::unordered_set< std::string > remove_set;
        //remove everything but the one forcing

        for(auto& itr : _metdata->stations())
        {
            if( itr->ID() != point_mode.forcing)
                remove_set.insert(itr->ID());
        }


        _metdata->prune_stations(remove_set);

        _outputs.erase(std::remove_if(_outputs.begin(),_outputs.end(),
                                      [this](const output_info& o){return o.name  != point_mode.output;}),
                       _outputs.end());


        LOG_DEBUG << _outputs.size();
        if ( _metdata->nstations() != 1 ||
                _outputs.size() !=1)
        {
            BOOST_THROW_EXCEPTION(model_init_error() << errstr_info(">1 station or outputs in point mode"));

        }
        for(auto s: _metdata->stations())
        {
            LOG_DEBUG << *s;
        }

        for(auto o:_outputs)
        {
            LOG_DEBUG << o.name;
        }
    }


    pt::json_parser::write_json((o_path  / "config.json" ).string(),cfg); // output a full dump of the cfg, after all modifications, to the output directory
    _cfg = cfg;

    LOG_DEBUG << "Finished initialization";

    LOG_DEBUG << "Determining module dependencies";
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
                BOOST_THROW_EXCEPTION(config_error() << errstr_info("Requested output " + var + " is not provided by any module."));
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



    LOG_DEBUG << "Allocating face variable storage";

    //we are going to make the assumption that every module can store face data.
    // However if this is onerous we can add a flag to the modules later
    std::set<std::string> module_list;
    for(auto& itr: _modules)
    {
        module_list.insert(itr.first->ID);
    }

    _mesh->init_face_data(_provided_var_module, _provided_var_vector, module_list);

    if(point_mode.enable)
    {
        for(auto itr:_chunked_modules)
        {
            for(auto jtr:itr)
            {
                if(jtr->parallel_type() == module_base::parallel::domain)
                {
                    BOOST_THROW_EXCEPTION(model_init_error() << errstr_info("Domain parallel module being run in point-mode."));
                }
            }
        }
    }


    timer c;
    LOG_DEBUG << "Running init() for each module";
    c.tic();


    for (auto& itr : _modules)
    {
        LOG_VERBOSE << itr.first->ID;
        itr.first->init(_mesh);
    }
    LOG_DEBUG << "Took " << c.toc<ms>() << "ms";

    //we do this here now because init is allowing a module to chance its mide and declar itself
    // data parallel or domain parallel after the fact.
    _schedule_modules();

//load a checkpoint as the last thing we do before a run
    if(_load_from_checkpoint  )
    {
        LOG_DEBUG << "Loading from checkpoint";
        c.tic();
        for (auto &itr : _chunked_modules)
        {
            //module calls
            for (auto &jtr : itr)
            {
                jtr->load_checkpoint(_mesh, _in_savestate);
            }
        }

        LOG_DEBUG << "Done loading snapshot [ " << c.toc<s>() << "s]";
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
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(
                        "An included json file cannot contain json sub-files."));
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
               BOOST_THROW_EXCEPTION(module_error() << errstr_info("Module " + m1->ID + " explicitly conflicts with " + m2->ID));
           }


            auto itr = std::find_first_of (m1_provides_var_names.begin(), m1_provides_var_names.end(),
                                           m2_provides_var_names.begin(), m2_provides_var_names.end());
            if( itr != m1_provides_var_names.end() )
            {
                overwrite_found=true;
                LOG_ERROR << "Module " << m1->ID << " and " << m2->ID << " both provide variable " << *itr;
            }
        }
    }

    if(overwrite_found)
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("A module's provides overwrites another module's provides. This is not allowed."));

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
            LOG_DEBUG << "Module [" << module->ID << "], no inter-module dependenices";
        } else
        {
            LOG_DEBUG << "Module [" << module->ID << "] checking against...";
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
	  auto itr_mod_provides_var_names = itr_module->get_variable_names_from_collection(*(itr_module->provides()));;
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
                        LOG_DEBUG << "\t\tAdding edge between " << module->ID << "[" << module->IDnum <<
                                  "] -> " << itr_module->ID << "[" << itr_module->IDnum << "] for var=" << *i <<
                                  std::endl;

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
                                LOG_WARNING << "Skipped adding edge between " << o.first << " and " << o.second <<
                                            " for var=" << *i << " because of user override" << std::endl;
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
                        LOG_DEBUG << "\t\tAdding optional edge between " << module->ID << "[" <<
                                  module->IDnum << "] -> " << itr_module->ID << "[" << itr_module->IDnum <<
                                  "] for var=" << *i << std::endl;

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
                                LOG_WARNING << "Skipped adding edge between " << o.first << " and " << o.second <<
                                            " for var=" << *i << " because of user override" << std::endl;
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
            LOG_ERROR << ss.str();
            BOOST_THROW_EXCEPTION(module_error() << errstr_info(ss.str()));
        }

        //check if our module has any met file dependenices
        if (module->depends_from_met()->size() == 0)
        {
            LOG_DEBUG << "Module [" << module->ID << "], no met file dependenices";
        } else
        {
            LOG_DEBUG << "Module [" << module->ID << "] has met file dependencies";
        }




        //check this modules met dependencies, bail if we are missing any.
        for (auto &depend_met_var : *(module->depends_from_met()))
        {
            auto i = std::find(_provided_var_met_files.begin(), _provided_var_met_files.end(), depend_met_var);
            if (i == _provided_var_met_files.end())
            {
                LOG_ERROR << "\t\t" << depend_met_var << "...[missing]";
                BOOST_THROW_EXCEPTION(module_error() << errstr_info("Missing dependency for " + depend_met_var));
            }
            LOG_DEBUG << "\t\t" << depend_met_var << "...[ok]";
        }

    }

    //great filter file for gvpr
    std::ofstream gvpr("filter.gvpr");;

    std::string font = "Helvetica";
    int fontsize = 11;

    std::string edge_str = "E[edgetype == \"%s\"] {\n color=\"/paired12/%i\";\n fontsize=%i;\n     fontname=\"%s\"\n }";
    int idx = 1;
    for (auto itr : graphviz_vars)
    {
        std::string edge = str_format(edge_str, itr.c_str(), idx, fontsize, font.c_str());
        idx++;
        gvpr << edge << std::endl;
    }

    gvpr.close();

    std::deque<int> topo_order;
    try
    {
        boost::topological_sort(g, std::front_inserter(topo_order));
    }
    catch(...)
    {

        std::ostringstream ssdot;
        boost::write_graphviz(ssdot, g, boost::make_label_writer(boost::get(&vertex::name, g)),
                              make_edge_writer(boost::get(&edge::variable, g)));
        std::string dot(ssdot.str());
        size_t pos = dot.find("\n");

        if (pos == std::string::npos)
            BOOST_THROW_EXCEPTION(config_error() << errstr_info("Unable to generate dot file"));
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
        ierr = std::remove("modules.dot"); CHK_SYSTEM_ERR(ierr);

        BOOST_THROW_EXCEPTION(config_error() << errstr_info("Module graph must be a DAG. Please review modules.pdf to determine where the cycle occured."));

    }




    std::ostringstream ssdot;
    boost::write_graphviz(ssdot, g, boost::make_label_writer(boost::get(&vertex::name, g)),
                          make_edge_writer(boost::get(&edge::variable, g)));
    std::string dot(ssdot.str());
    size_t pos = dot.find("\n");

    if (pos == std::string::npos)
        BOOST_THROW_EXCEPTION(config_error() << errstr_info("Unable to generate dot file"));

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
    int ierr = std::system("gvpr -c -f filter.gvpr -o modules.dot modules.dot.tmp"); CHK_SYSTEM_ERR(ierr);
    ierr = std::system("dot -Tpdf modules.dot -o modules.pdf"); CHK_SYSTEM_ERR(ierr);
    ierr = std::remove("modules.dot.tmp"); CHK_SYSTEM_ERR(ierr);
    ierr = std::remove("filter.gvpr"); CHK_SYSTEM_ERR(ierr);

    std::stringstream ss;
    size_t order = 0;
    for (std::deque<int>::const_iterator i = topo_order.begin(); i != topo_order.end(); ++i)
    {
        ss << _modules.at(*i).first->ID << "->";
        _modules.at(*i).second = order;
        order++;
    }

    std::string s = ss.str();
    LOG_DEBUG << "Build order: " << s.substr(0, s.length() - 2);


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
    LOG_DEBUG << "_modules order after sort: " << s.substr(0, s.length() - 2);

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
        LOG_DEBUG << "Chunking module: " << itr.first->ID;
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
        LOG_DEBUG << "Chunk " << (itr.at(0)->parallel_type() == module_base::parallel::data ? "data" : "domain") <<
                  " " << chunks << ": ";
        for (auto &jtr : itr)
        {
            LOG_DEBUG << jtr->ID;
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


    LOG_DEBUG << "Loading first timestep's met data";
    // Populate the stations with the first timestep's data.
    // We can do this _once_ without incrementing the internal iterators
    _metdata->next();

    LOG_DEBUG << "Starting model run";

    c.tic();

    double meantime = 0;
    size_t current_ts = 0;
    _global->timestep_counter = 0; //use this to pass the timestep info to the modules for easier debugging specific timesteps
    size_t max_ts = _metdata->n_timestep();
    bool done = false;

        while (!done)
        {
            boost::posix_time::ptime t;

            _global->_current_date = _metdata->current_time();

            LOG_DEBUG << "Timestep: " << _global->posix_time() << "\tstep#"<<current_ts;

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
                        for (size_t i = 0; i < _mesh->size_faces(); i++)
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
                LOG_ERROR << "Exception at timestep: " << _global->posix_time();
                //if we die in a module, try to dump our time series out so we can figure out wtf went wrong
                LOG_ERROR << "Exception has occured. Timeseries and meshes WILL BE INCOMPLETE!";
                *_end_ts = _global->posix_time();
                done = true;
                LOG_ERROR << boost::diagnostic_information(e);

            }
            catch(std::exception& e)
            {
                LOG_ERROR << "Exception at timestep: " << _global->posix_time();
                LOG_ERROR << "Unknown exception:";
                LOG_ERROR << e.what();
                *_end_ts = _global->posix_time();
                done = true;
            }

            //check that we actually need a mesh output.
            for (auto &itr : _outputs)
            {
                if(itr.type == output_info::output_type::mesh)
                {
                    std::vector<std::string> output;
                    output.assign(itr.variables.begin(),itr.variables.end()); //convert to list to match internal lists

                    _mesh->update_vtk_data(output); //update the internal vtk mesh
                    break; // we're done as soon as we've called update once. No need to do it multiple times.
                }
            }

            // save the current state
            if(_do_checkpoint && (current_ts % _checkpoint_feq ==0) )
            {
                LOG_DEBUG << "Checkpointing...";
                c.tic();
                for (auto &itr : _chunked_modules)
                {
                    //module calls
                    for (auto &jtr : itr)
                    {
                        jtr->checkpoint(_mesh, _savestate);
                    }
                }
                std::stringstream timestr;

                timestr << _global->posix_time() + boost::posix_time::seconds(_global->_dt); // start from current TS + dt

                //also write it out in seconds because netcdf is struggling with the string
                unsigned long long int ts_sec = _global->posix_time_int()+_global->_dt;

                _savestate.get_ncfile().putAtt("restart_time",timestr.str());
                _savestate.get_ncfile().putAtt("restart_time_sec", netCDF::ncUint64,ts_sec);

                LOG_DEBUG << "Done checkpoint [ " << c.toc<s>() << "s]";
            }

            for (auto &itr : _outputs)
            {
                if (itr.type == output_info::output_type::mesh)
                {
                    // check if we should output or not
                    bool should_output = false;

                    if(itr.only_last_n != -1)
                    {
                        auto ts_left = max_ts - current_ts;
                        if( ts_left <= itr.only_last_n) // if we are within the last n timesteps, output
                            should_output = true;
                    }
                    else
                    {
                        if(current_ts % itr.frequency == 0)
                            should_output = true;
                    }



                    if(should_output)
                    {

                        #pragma omp parallel
                        {
                            #pragma omp single
                            {
                                for (auto jtr : itr.mesh_output_formats)
                                {
                                    #pragma omp task
                                    {
                                        std::string base_name = itr.fname + std::to_string(_global->posix_time_int());
                                        boost::filesystem::path p(base_name);

                                        if (jtr == output_info::mesh_outputs::vtu  )
                                        {

                                            // this really only works if we let rank0 handle the io.
                                            // If we let each process do it, they walk all over each other's output
#ifdef USE_MPI
                                            if(_comm_world.rank() == 0)
                                            {
                                                for(int rank = 0; rank < _comm_world.size(); rank++)
                                                {
#else
                                                    int rank = 0;
#endif
                                                    pt::ptree &dataset = pvd.add("VTKFile.Collection.DataSet", "");
                                                    dataset.add("<xmlattr>.timestep", _global->posix_time_int());
                                                    dataset.add("<xmlattr>.group", "");
                                                    dataset.add("<xmlattr>.part", rank);
                                                    dataset.add("<xmlattr>.file", p.filename().string()+"_"+std::to_string(rank) + ".vtu");
#ifdef USE_MPI
                                                }
                                            }
#endif

                                            //because a full path can be provided for the base_name, we need to strip this off
                                            //to make it a relative path in the xml file.

#ifdef USE_MPI
                                            _mesh->write_vtu(base_name + "_"+std::to_string(_comm_world.rank() )+ ".vtu");
#else
                                            _mesh->write_vtu(base_name + "_"+std::to_string(rank)+ ".vtu");
#endif

                                        }
                                    }
                                }
                            }
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
                    for (auto v : _provided_var_module)
                    {
                        auto data = (*itr.face)[v];
                        itr.ts.at(v, current_ts) = data;
                    }
                }
            }

            if(!_metdata->next())
                done = true;

            auto timestep = c.toc<ms>();
            meantime += timestep;

            current_ts++;
            _global->timestep_counter++;

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

            LOG_DEBUG << "Avg timestep duration " << s << "\tEstimated completion: " << boost::posix_time::to_simple_string(pt);
            _global->first_time_step = false;


        }
        double elapsed = c.toc<s>();
        LOG_DEBUG << "Total runtime was " << elapsed << "s";



    std::string base_name="";

    for (auto &itr : _outputs)
    {
        if (itr.type == output_info::output_type::mesh)
        {

#ifdef USE_MPI
            if(_comm_world.rank() == 0)
            {
#endif
#if (BOOST_VERSION / 100 % 1000) < 56
                pt::write_xml(base_name + ".pvd",
                              pvd, std::locale(), pt::xml_writer_make_settings<char>(' ', 4));
#else
                pt::write_xml(itr.fname + ".pvd",
                              pvd, std::locale(), pt::xml_writer_settings<std::string>(' ', 4));
                break;
#endif
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
    }


    if(_notification_script != "")
    {
        LOG_DEBUG << "Calling notification script";
        int ierr = std::system(_notification_script.c_str()); CHK_SYSTEM_ERR(ierr);
    }
}

void core::end(const bool abort)
{
#ifdef USE_MPI
    if(abort)
    {
        LOG_ERROR << "An exception has occurred, requesting MPI Abort!";
        _mpi_env.abort(-1);
    }
#endif

    LOG_DEBUG << "Cleaning up";
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

    LOG_DEBUG << "Populating each face's station list";

    for (size_t i = 0; i < _mesh->size_faces(); i++)
    {
        auto f = _mesh->face(i);

        if ( f->stations().size() == 0 )
        {
            auto stations = _metdata->get_stations(f->get_x(), f->get_y());
            f->stations().insert(std::end(f->stations()), std::begin(stations), std::end(stations));

            f->nearest_station() = _metdata->nearest_station(f->get_x(), f->get_y()).at(0);

        }  else
        {
            CHM_THROW_EXCEPTION(mesh_error,"Face station list already populated.");
        }
    }

}

void core::populate_distributed_station_lists()
{
    std::vector<std::shared_ptr<station>> _stations; //stations to cull

    using th_safe_multicontainer_type = std::vector< std::shared_ptr<station> >[];
    std::unique_ptr< th_safe_multicontainer_type > th_local_stations;
    std::vector< std::shared_ptr<station> > mpi_local_stations;

    LOG_DEBUG << "Populating each MPI process's station list";

#pragma omp parallel
    {
        // We want an array of vectors, so that OMP threads can increment them
        // separately, then join them afterwards
#pragma omp single
        {
            th_local_stations = std::make_unique< th_safe_multicontainer_type >(omp_get_num_threads());
        }
#pragma omp for
        for(size_t face_index=0; face_index< _mesh->size_faces(); ++face_index)
        {
            // face_index is a local index... get the face handle
            auto face = _mesh->face(face_index);
            if ( face->stations().empty() )
            { // only perform if faces' stationlists are set
                BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Face station lists must be populated before populating distributed MPI station lists."));
            }
            for (auto &p : face->stations())
            {
                th_local_stations[omp_get_thread_num()].push_back(p);
            }
        }
        // Join the vectors via a single thread in t operations
        //  NOTE future optimizations:
        //   - reserve space for insertions into mpi_local_stations
        //   - can be done recursively in log2(num_threads) operations
#pragma omp single
        {
            for(int thread_idx=0;thread_idx<omp_get_num_threads();++thread_idx)
            {
                mpi_local_stations.insert(std::end(mpi_local_stations),
                                          std::begin(th_local_stations[thread_idx]), std::end(th_local_stations[thread_idx]));
            }
        }

    }

    // Remove duplicates by converting to a set
    std::unordered_set< std::shared_ptr<station>  > keep_set(std::begin(mpi_local_stations), std::end(mpi_local_stations));

    std::unordered_set< std::string > remove_set;
    for(auto& itr: _metdata->stations())
    {
        if( keep_set.find(itr) == keep_set.end() ) // not found in the set we want to keep, mark for removal
            remove_set.insert(itr->ID());
    }

    // Store the local stations in the triangulations mpi-local stationslist vector
    _metdata->prune_stations(remove_set);

#ifdef USE_MPI
    LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _metdata->nstations() << " locally owned stations.";
#endif
}
