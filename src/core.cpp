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

#include "core.h"

core::core()
{
    BOOST_LOG_FUNCTION();

    _enable_ui=false; //default to no UI
    _start_ts = nullptr;
    _end_ts = nullptr;
    _per_triangle_timeseries = false;
    _interpolation_method = interp_alg::tpspline;

    //default logging level
    _log_level = debug;
    _output_station_ptv = true;
    _use_netcdf=false;
    _load_from_checkpoint=false;
    _do_checkpoint=false;
}

core::~core()
{
    LOG_DEBUG << "Finished";
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

    //enable/disable ncurses UI. default is enable
    boost::optional<bool> u = value.get_optional<bool>("ui");
    if (u)
    {
        _enable_ui = *u;
        LOG_DEBUG << "Set ui to " << *u;
    }

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


    // project name
    boost::optional<std::string> prj = value.get_optional<std::string>("prj_name");
    if (prj)
    {
        //_prj_name = *prj;
        _ui.write_model_name(*prj);
        LOG_DEBUG << "Set project name to " << *prj;
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

    // full time series per triangle.
    boost::optional<bool> per_triangle_storage = value.get_optional<bool>("per_triangle_timeseries");
    if (per_triangle_storage)
    {
        _per_triangle_timeseries = *per_triangle_storage;
        LOG_DEBUG << "per triangle timer series storage: " << _per_triangle_timeseries;
    }
    else
    {
        _per_triangle_timeseries = false;
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
        _global->station_search_radius = *radius;
        _global->get_stations = boost::bind( &global::get_stations_in_radius,_global,_1,_2, *radius);
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

        _global->N = n;
        _global->get_stations = boost::bind( &global::nearest_station,_global,_1,_2, n);
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
            //    LOG_DEBUG << "inserted module " << module;
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

        std::string module = itr;
        LOG_DEBUG << "Module ID=" << module;


        //try grabbing a config for this module, empty string default
        pt::ptree cfg;
        try
        {
            cfg = config.get_child(module);
        } catch (pt::ptree_bad_path &e)
        {
            LOG_DEBUG << "No config for " << module;
        }
        boost::shared_ptr<module_base> m(_mfactory.get(module, cfg));
        //internal tracking of module initialization order
        m->IDnum = modnum;

        //assign the internal global param pointer to our global
        m->global_param = _global;

        modnum++;
        _modules.push_back(
                std::make_pair(m, 1)); //default to 1 for make ordering, we will set it later in determine_module_dep
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
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkStringArray> labels = vtkSmartPointer<vtkStringArray>::New();
    labels->SetName("Station name");
    std::vector< std::pair<std::string, pt::ptree> > forcings;

    //positive offset going west. So the normal UTC-6 would be UTC_offset:6
    _global->_utc_offset = value.get("UTC_offset",0);
    LOG_DEBUG << "Applying UTC offset to ALL forcing files. UTC+" << std::to_string(_global->_utc_offset);

    _find_and_insert_subjson(value);

    //need to determine if we have been given a netcdf file
    _use_netcdf = value.get("use_netcdf",false);

    size_t nstations = 0;
    timer c;

    //we need to treat this very differently than the txt files
    if(_use_netcdf)
    {
        LOG_DEBUG << "Found NetCDF file";
        std::string file = value.get<std::string>("file");

        // core needs to setup all the stations as timeseries doens't know what to do with a full netcdf file.

        try
        {
            nc.open_GEM(file);
        }catch(netCDF::exceptions::NcException& e)
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(e.what()));
        }


        nstations = nc.get_xsize() * nc.get_ysize();
    } else
    {
        for (auto &itr : value)
        {
            if(itr.first != "UTC_offset")
            {
                std::string station_name = itr.first.data();
                forcings.push_back(make_pair(station_name, value.get_child(station_name)));
                ++nstations;
            }
        }
        LOG_DEBUG << "Found # stations = " <<  nstations;

        if(nstations == 0)
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("No input forcing files found!"));
        }
    }


    _global->_stations.resize(nstations);
    tbb::concurrent_vector<  boost::shared_ptr<station> > pstations;
    pstations.resize(nstations);

    OGRSpatialReference insrs;

    if(_use_netcdf)
    {
        insrs.SetWellKnownGeogCS("WGS84");
    } else
    {
        insrs.SetWellKnownGeogCS("WGS84");
    }


    OGRSpatialReference outsrs;
    auto err = outsrs.importFromProj4(_mesh->proj4().c_str());

    OGRCoordinateTransformation* coordTrans =  nullptr;

    if(!_mesh->is_geographic())
        coordTrans = OGRCreateCoordinateTransformation(&insrs, &outsrs);

    labels->SetNumberOfValues(nstations);

    LOG_DEBUG << "Reading forcing file data";
    c.tic();

    if(_use_netcdf)
    {
        try
        {
            auto variables = nc.get_variable_names();
            auto date_vec = nc.get_datevec();


            LOG_DEBUG << "Grid is (y)" << nc.get_ysize() << " by (x)" << nc.get_xsize();


            LOG_DEBUG << "Loading lat/long grid...";
            c.tic();

            auto lat = nc.get_lat();
            auto lon = nc.get_lon();
            auto e = nc.get_z();
            LOG_DEBUG << "Took " << c.toc<s>() << "s";

            LOG_DEBUG << "Initializing datastructure";
            c.tic();

// #pragma omp parallel for
// hangs, unclear why, critical sections around the json and gdal calls don't seem to help
            for (size_t y = 0; y < nc.get_ysize(); y++)
            {
                for (size_t x = 0; x < nc.get_xsize(); x++)
                {


                    double latitude = 0;
                    double longitude = 0 ;
                    double z = 0;

                    latitude = lat[y][x];;
                    longitude = lon[y][x];
                    z = e[y][x];

//                    latitude = nc.get_lat(x, y);
//                    longitude = nc.get_lon(x, y);
//                    z = nc.get_z(x, y);


//                    LOG_DEBUG << "lat=" << latitude << "\tlong="<<longitude<<"\tz="<<z;
//                    LOG_DEBUG << "lat=" << lat[y][x] << "\tlong="<<lon[y][x]<<"\tz="<<e[y][x];
//                    LOG_DEBUG << "------------------------";


                    size_t index = x + y * nc.get_xsize();
                    std::string station_name = std::to_string(index); // these don't really have names

                    //project mesh, need to convert the input lat/long into the coordinate system our mesh is in
                    if (!_mesh->is_geographic())
                    {
                        if (!coordTrans->Transform(1, &longitude, &latitude))
                        {
                            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(
                                    "Station=" + station_name + ": unable to convert coordinates to mesh format."));
                        }
                    }

                    double elevation = z;

                    //this ctor will init an empty timeseries for us
                    boost::shared_ptr<station> s = boost::make_shared<station>(station_name, longitude, latitude,
                                                                               elevation);
                    s->raw_timeseries()->init(variables, date_vec);
                    s->reset_itrs();

                    try
                    {
                        auto filter_section = value.get_child("filter");

                        for (auto &jtr: filter_section)
                        {
                            auto name = jtr.first.data();

                            boost::shared_ptr<filter_base> filter(_filtfactory.get(name, jtr.second));
                            filter->init(s);
                            _netcdf_filters[name] = filter;

                        }

                    } catch (pt::ptree_bad_path &e)
                    {
                        //ignore bad path, means we don't have a filter
                    }

                    //index this linear array as if it were 2D to make the lazy load in the main run() loop easier.
                    //it will allow us to pull out the station for a specific x,y more easily.
                    pstations.at(index) = s;

                }
            }
            LOG_DEBUG << "Took " << c.toc<s>() << "s";

        } catch(netCDF::exceptions::NcException& e)
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info(e.what()));
        }

    } else
    {
        if (nstations != forcings.size())
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Something has gone wrong in forcing file read."));
        }

        //#pragma omp parallel for
        //TODO: this dead locks, not sure why <-- json reader isn't parallel safe
        for(size_t i =0; i < nstations; ++i)
        {
            auto& itr = forcings.at(i);

            std::string station_name = itr.first;//.data(); n

            boost::shared_ptr<station> s = boost::make_shared<station>();

            s->ID(station_name);

            double longitude=0;
            double latitude=0;
            try
            {
                longitude = itr.second.get<double>("longitude");
                latitude = itr.second.get<double>("latitude");
            }catch(boost::property_tree::ptree_bad_path& e)
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Station " + station_name + " missing lat/long coordinate."));

            }

            if( (latitude > 90 || latitude < -90) ||
                (longitude > 180 || longitude < -180) )
            {
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Station " + station_name + " coordinate is invalid."));
            }

            //project mesh, need to convert the input lat/long into the coordinate system our mesh is in
            if(!_mesh->is_geographic())
            {
                if(!coordTrans->Transform(1, &longitude, &latitude))
                {
                    BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Station=" + station_name + ": unable to convert coordinates to mesh format."));
                }
            }

            s->x(longitude);
            s->y(latitude);

            double elevation = itr.second.get<double>("elevation");
            s->z(elevation);

            std::string file = itr.second.get<std::string>("file");
            auto f = cwd_dir / file;
            s->open(f.string());

            //filters with text timeseries behave differently than filters with the netcdf
            //these will be run once, over the entire timeseries. netcdf will call the filters and run it on every lazyload.
            try
            {
                auto filter_section = itr.second.get_child("filter");

                for (auto &jtr: filter_section)
                {
                    auto name = jtr.first.data();

                    boost::shared_ptr<filter_base> filter(_filtfactory.get(name, jtr.second));
                    filter->init(s);

                    //save this filter to run later
                    //each station is associated with a list of filters to run
                    _txtmet_filters[s->ID()].push_back(filter);

                }

            } catch (pt::ptree_bad_path &e)
            {
                //ignore
            }

            pstations.at(i) = s;

        }

    }

    for(size_t i = 0; i<nstations;i++)
    {
        auto s = pstations.at(i);
        if(_output_station_ptv)
        {
            points->InsertNextPoint(s->x(), s->y(), s->z());
            labels->SetValue(i, s->ID() );
        }
        auto cf = _mesh->find_closest_face(s->x(),s->y());
        s->set_closest_face(cf->cell_id);

        //do a few things behind _global's back for efficiency.
//        auto s = pstations.at(i);
        _global->_stations.at(i) = s;
        _global->_dD_tree.insert(boost::make_tuple(global::Kernel::Point_2(s->x(), s->y()), s));
    }


    delete coordTrans;

    LOG_DEBUG << "Took " << c.toc<s>() << "s";

    LOG_DEBUG << "Finished reading stations";

    if(_output_station_ptv)
    {

        vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        polydata->GetPointData()->AddArray(labels);
        // Write the file
        vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        auto f = o_path / "stations.vtp";
        writer->SetFileName(f.string().c_str());

#if VTK_MAJOR_VERSION <= 5
        writer->SetInput(polydata);
#else
        writer->SetInputData(polydata);
#endif

        writer->Write();
    }


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

    pt::ptree mesh = read_json(mesh_path);

    //we need to check if we've read in a geographic (lat/long) mesh or a UTM mesh. We then need to swap in the right distance and point_bearing functions
    //so the modules and future code can blindly use them without worrying about these things

    bool is_geographic = (bool)mesh.get<int>("mesh.is_geographic");
    _global->_is_geographic = is_geographic; // save it here so modules can determine if this is true
    if(is_geographic)
    {

        math::gis::point_from_bearing = & math::gis::point_from_bearing_latlong;
        math::gis::distance = &math::gis::distance_latlong;
    } else
    {

        math::gis::point_from_bearing = &math::gis::point_from_bearing_UTM;
        math::gis::distance = &math::gis::distance_UTM;
    }

    bool triarea_found = false;
    //see if we have additional parameter files to load
    try
    {
        for(auto &itr : value.get_child("parameters"))
        {
            LOG_DEBUG << "Parameter file: " << itr.second.data();

            auto param_mesh_path = (cwd_dir / itr.second.data()).string();

            pt::ptree param_json = read_json(param_mesh_path);

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

    }
    catch(pt::ptree_bad_path &e)
    {
        LOG_DEBUG << "No addtional parameters found in mesh section.";
    }
    if(is_geographic && !triarea_found)
    {
        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Geographic meshes require the triangle area be present in a .param file. Please include this."));
    }

    try
    {
        for(auto &itr : value.get_child("initial_conditions"))
        {
            LOG_DEBUG << "Initial condition file: " << itr.second.data();

            auto param_mesh_path = (cwd_dir / itr.second.data()).string();

            pt::ptree ic_json = read_json(param_mesh_path);

            for(auto& ktr : ic_json)
            {
                //use put to ensure there are no duplciate parameters...
                std::string key = ktr.first.data();
                mesh.put_child( "initial_conditions." + key ,ktr.second);
                LOG_DEBUG << "Inserted initial condition " << ktr.first.data() << " into the config tree.";
            }

        }
    }
    catch(pt::ptree_bad_path &e)
    {
        LOG_DEBUG << "No addtional initial conditions found in mesh section.";
    }

    _provided_parameters = _mesh->from_json(mesh);

    if (_mesh->size_faces() == 0)
        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Mesh size = 0!"));

    _ui.write_mesh_details(_mesh->size_faces());

    LOG_DEBUG << "Initializing DEM mesh attributes";
#pragma omp parallel for
    for (size_t i = 0; i < _mesh->size_faces(); i++)
    {
        auto face = _mesh->face(i);
        face->slope();
        face->aspect();
        face->center();
        face->normal();
    }



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
                BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Missing latitude and/or longitude for " + out.name));
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

                if(!coordTrans->Transform(1, &out.longitude, &out.latitude))
                {
                    BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Output=" + out.name + ": unable to convert coordinates to mesh format."));
                }

                delete coordTrans;
            }



            out.face = _mesh->locate_face(out.longitude, out.latitude);


            if(out.face == nullptr)
                BOOST_THROW_EXCEPTION(config_error() <<
                                                     errstr_info(
                                                             "Requested an output point that is not in the triangulation domain. Pt:"
                                                             + std::to_string(out.longitude) + "," + std::to_string(out.latitude) + " name: " + out.name));

            //set the point to be the center of the triangle that the output point lies on
            points->InsertNextPoint ( out.face->get_x(), out.face->get_y(), out.face->get_z() );
            labels->InsertNextValue(out.name);

            LOG_DEBUG << "Triangle geometry for output triangle = " << out_type << " slope: " << out.face->slope() * 180./M_PI << " aspect:" << out.face->aspect() * 180./M_PI;

            out.face->_debug_name = out.name; //out_type holds the station name
            out.face->_debug_ID = ID;
            ++ID;
        }
        else if (out_type == "mesh")
        {
            out.type = output_info::mesh;

            auto fname = itr.second.get<std::string>("base_name","output");
            auto f = msh_path / fname;
            boost::filesystem::create_directories(f.parent_path());
            out.fname = f.string();

            _mesh->write_param_to_vtu( itr.second.get("write_parameters",true) ) ;

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

            out.mesh_output_formats.push_back(output_info::mesh_outputs::vtu);


            // can do all non-vtu formats with vtu2geo tool, which is also faster. So move towards removing all this
            //functionality

//            for (auto &jtr: itr.second.get_child("format"))
//            {
//                LOG_DEBUG << "Output format found: " << jtr.second.data();
//                if (jtr.second.data() == "vtu")
//                    out.mesh_output_formats.push_back(output_info::mesh_outputs::vtu);
//                if (jtr.second.data() == "vtp")
//                    out.mesh_output_formats.push_back(output_info::mesh_outputs::vtp);
//            }

        } else
        {
            LOG_WARNING << "Unknown output type: " << itr.second.data();
        }

        _outputs.push_back(out);
    }
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(labels);
    // Write the file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    //output this to the same folder as the points are written out to
    auto f = pts_path / "output_points.vtp";
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
    std::string version = "CHM 0.1 " GIT_BRANCH "/" GIT_COMMIT_HASH;

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
        exit(1);
    }
    else if (vm.count("version"))
    {
        cout << version << std::endl;
        exit(1);
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
    cwd_dir = path.parent_path();

    std::string log_dir = "log";
    auto log_path = cwd_dir / log_dir;
    std::string log_name = "CHM_" + log_start_time + ".log";

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
                    expr::format("%1% %2% [%3%]: %4%")
                    % expr::attr<boost::posix_time::ptime>("TimeStamp")
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

    logging::core::get()->add_global_attribute("TimeStamp", attrs::local_clock());
    logging::core::get()->add_global_attribute("Scope", attrs::named_scope());



    logging::core::get()->add_sink(_log_sink);

    LOG_DEBUG << "Logger initialized. Writing to cout and " + log_name;


#ifdef NOMATLAB
    _engine = boost::make_shared<maw::matlab_engine>();


    _engine->start();
    _engine->set_working_dir();

    LOG_DEBUG << "Matlab engine started";
#endif

    _global = boost::make_shared<global>();

    //don't just abort and die
    gsl_set_error_handler_off();


#ifdef _OPENMP
    LOG_DEBUG << "Built with OpenMP support, #threads = " << omp_get_max_threads();
#endif


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
     *  modules
     *  meshes
     *  forcing
     * The rest may be optional, and will override the defaults.
     */
    config_modules(cfg.get_child("modules"), cfg.get_child("config"), cmdl_options.get<3>(), cmdl_options.get<4>());
    config_meshes(cfg.get_child("meshes")); // this must come before forcing, as meshes initializes the required distance functions based on geographic/utm meshes

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


    //if we run under a shitty terminal that doesn't support ncurses, or GDB
    //we do have to turn this off and fall back to just showing cout
    if(_enable_ui)
    {
        try
        {
            _ui.init();
        } catch (model_init_error &e)
        {
            LOG_WARNING << "ncurses did not init, falling back to cout";
            _enable_ui = false;
        }
    }


    boost::filesystem::path full_path(boost::filesystem::current_path());
    _ui.write_cwd(full_path.string());



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

        //remove everything but the one forcing
        _global->_stations.erase(std::remove_if(_global->_stations.begin(),_global->_stations.end(),
                       [this](boost::shared_ptr<station> s){return s->ID() != point_mode.forcing;}),
                                _global->_stations.end());

        _outputs.erase(std::remove_if(_outputs.begin(),_outputs.end(),
                                      [this](output_info o){return o.name  != point_mode.output;}),
                       _outputs.end());

        LOG_DEBUG << _global->_stations.size();
        LOG_DEBUG << _outputs.size();
        if ( _global->_stations.size() != 1 ||
                _outputs.size() !=1)
        {
            BOOST_THROW_EXCEPTION(model_init_error() << errstr_info(">1 station or outputs in point mode"));

        }
        for(auto s: _global->_stations)
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

    LOG_DEBUG << "Init variables mapping";
    _global->_variables.init_from_file("Un-init path");

    LOG_DEBUG << "Determining module dependencies";
    _determine_module_dep();

    LOG_DEBUG << "Initializing and allocating memory for timeseries";

    if (_global->_stations.size() == 0)
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("no stations"));


    auto start_time = _global->_stations.at(0)->date_timeseries().at(0);
    auto end_time = _global->_stations.at(0)->date_timeseries().back();

    //if the user gives both a custom start and end time, don't bother with this.
    //this does assume the user knows what they are doing and subsets a time period that is correct
    //if they do not, the sanity check after this period that double checks consistent timestepping will catch it
    if(!_start_ts && !_end_ts)
    {
        // find the latest start time and the earliest end time
        for (size_t i = 0; i < _global->_stations.size(); ++i)
        {
            auto tmp = _global->_stations.at(i)->date_timeseries().at(0);

            if (tmp > start_time)
                start_time = tmp;

            tmp = _global->_stations.at(i)->date_timeseries().back();
            if (tmp < end_time)
                end_time = tmp;
        }
    }

    // If we ended on time T, restart from T+1. T+1 is written out to attr, so we can just start from this
    if(_load_from_checkpoint)
    {
        //better but is loading incorrectly, seems to be an off by one?
//        std::string timestr;
//        _in_savestate.get_ncfile().getAtt("restart_time").getValues(&timestr);

        size_t t = 0;
        _in_savestate.get_ncfile().getAtt("restart_time_sec").getValues(&t);

///        _start_ts = new boost::posix_time::ptime(boost::posix_time::time_from_string(timestr));
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

    if( _global->_stations.at(0)->date_timeseries().size() == 1)
    {
        BOOST_THROW_EXCEPTION(model_init_error() << errstr_info("Unable to determine model timestep from only 1 input timestep."));
    }
    //figure out what our timestepping is. Needs to happen before we subset as we may end up with only 1 timestep
    auto t0 = _global->_stations.at(0)->date_timeseries().at(0);
    auto t1 = _global->_stations.at(0)->date_timeseries().at(1);
    auto dt = (t1 - t0);
    _global->_dt = dt.total_seconds();
    LOG_DEBUG << "model dt = " << _global->dt() << " (s)";


    LOG_DEBUG << "Subsetting station timeseries";

#pragma omp for
    for(size_t i = 0; i < _global->_stations.size(); ++i)
    {
        _global->_stations.at(i)->raw_timeseries()->subset(*_start_ts, *_end_ts);
        _global->_stations.at(i)->reset_itrs();
        auto s= _global->_stations.at(i);
         LOG_VERBOSE << s->ID() << " Start = " << s->date_timeseries().front() << " End = " << s->date_timeseries().back();
    }

    //ensure all the stations have the same start and end times
    // per-timestep agreeent happens during runtime.
    start_time = _global->_stations.at(0)->date_timeseries().at(0);
    end_time = _global->_stations.at(0)->date_timeseries().back();


    for (size_t i = 1; //on purpose to skip first station
         i < _global->_stations.size();
         i++)
    {
        if (_global->_stations.at(i)->date_timeseries().at(0) != start_time ||
            _global->_stations.at(i)->date_timeseries().back() != end_time)
        {
            BOOST_THROW_EXCEPTION(forcing_timestep_mismatch()
                                  <<
                                  errstr_info("Timestep mismatch at station: " + _global->_stations.at(i)->ID()));
        }
    }

    //set interpolation algorithm
    _global->interp_algorithm = _interpolation_method;// interp_alg::tpspline;




    //setup output timeseries sinks
    auto date = _global->_stations.at(0)->date_timeseries();

    for (auto &itr : _outputs)
    {
        if (itr.type == output_info::output_type::time_series)
        {
            itr.ts.init(_provided_var_module, date); /*length of all the vectors to initialize*/
        }
    }



    #pragma omp parallel for
    for (size_t it = 0; it < _mesh->size_faces(); it++)
    {
        Delaunay::Face_handle face = _mesh->face(it);
        if(point_mode.enable && face->_debug_name != _outputs[0].name )
            continue;

        //_per_triangle_timeseries=true will create a full timeseries on each triangle. this takes up a crazy amount of ram
        if (_per_triangle_timeseries)
        {
            face->init_time_series(_provided_var_module, date); /*length of all the vectors to initialize*/
        }

        if (!_per_triangle_timeseries)
        {
            //only do 1 init. The time is wrong, but that is ok as it's never used
            timeseries::date_vec d;
            d.push_back(date[0]);
            face->init_time_series(_provided_var_module, d);
        }

    }

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

            ompException oe;
            oe.Run([&]
                   {
                       itr.first->init(_mesh);
                   });
            oe.Rethrow();
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

    std::list<std::string> all_provides;

    bool overwrite_found = false;
    //make sure modules cannot override another module's output
    for (auto &m1 : _modules)
    {
        for (auto &m2 : _modules)
        {
            if(m1.first->ID == m2.first->ID )
                continue;

            auto itr = std::find_first_of (m1.first->provides()->begin(), m1.first->provides()->end(),
                                           m2.first->provides()->begin(), m2.first->provides()->end());
            if( itr != m1.first->provides()->end() )
            {
                overwrite_found=true;
                LOG_ERROR << "Module " << m1.first->ID << " and " << m2.first->ID << " both provide variable " << *itr;
            }
        }
    }

    if(overwrite_found)
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("A module's provides overwrites another module's provides. This is not allowed."));

    //build a list of variables provided by the met files, culling duplicate variables from multiple stations.
    for (size_t i = 0; i < _global->_stations.size(); i++)
    {
        auto vars = _global->_stations.at(i)->list_variables();
        _provided_var_met_files.insert(vars.begin(), vars.end());
    }

    //loop through each module
    for (auto &module : _modules)
    {
        //Generate a  list of all variables,, provided from this module, and append to the total list, culling duplicates between modules.
        _provided_var_module.insert(module.first->provides()->begin(), module.first->provides()->end());

//        vertex v;
//        v.name = module.first->ID;
//        boost::add_vertex(v,g);
        g[module.first->IDnum].name = module.first->ID;
//        output_graph << module.first->IDnum << " [label=\"" << module.first->ID  << "\"];" << std::endl;

        //names.at(module.first->IDnum) =  module.first->ID;

        //check intermodule depends
        if (module.first->depends()->size() == 0)  //check if this module requires dependencies
        {
            LOG_DEBUG << "Module [" << module.first->ID << "], no inter-module dependenices";
        } else
        {
            LOG_DEBUG << "Module [" << module.first->ID << "] checking against...";
        }


        //make a copy of the depends for this module. we will then count references to each dependency
        //this will allow us to know what is missing
        std::map<std::string, size_t> curr_mod_depends;

        //populate a list of all the dependencies
        for (auto &itr : *(module.first->depends()))
        {
            curr_mod_depends[itr] = 0; //ref count init to 0
        }

        //iterate over all the modules,
        for (auto &itr : _modules)
        {
            //don't check against our module
            if (module.first->ID.compare(itr.first->ID) != 0)
            {
                //loop through each required variable of our current module
                for (auto &depend_var : *(module.first->depends()))
                {
                    //LOG_DEBUG << "\t\t[" << itr.first->ID << "] looking for var=" << depend_var;

                    auto i = std::find(itr.first->provides()->begin(), itr.first->provides()->end(), depend_var);
                    if (i != itr.first->provides()->end()) //itr provides the variable we are looking for
                    {
                        LOG_DEBUG << "\t\tAdding edge between " << module.first->ID << "[" << module.first->IDnum <<
                                  "] -> " << itr.first->ID << "[" << itr.first->IDnum << "] for var=" << *i <<
                                  std::endl;

                        //add the dependency from module -> itr, such that itr will come before module
                        edge e;
                        e.variable = *i;

                        //check if we need to ignore this edge as a result of user specific overrdige
                        //this will avoid a cycle is this exists.
                        bool ignore = false;
                        for (auto o : _overrides)
                        {
                            if (o.first == module.first->ID &&
                                o.second == itr.first->ID)
                            {
                                ignore = true;
                                LOG_WARNING << "Skipped adding edge between " << o.first << " and " << o.second <<
                                            " for var=" << *i << " because of user override" << std::endl;
                            }
                        }

                        if (!ignore)
                            boost::add_edge(itr.first->IDnum, module.first->IDnum, e, g);

                        //even if we ignore, inc our depencies so we don't fail later

                        //output_graph << itr.first->IDnum << "->" << module.first->IDnum << " [label=\"" << *i << "\"];" << std::endl;
                        curr_mod_depends[*i]++; //ref count our variable

                        //only link in the output if we aren't ignoring
                        if (!ignore)
                            graphviz_vars.insert(*i);
                    }
                }

                //loop through each required variable of our current module
                for (auto &optional_var : *(module.first->optionals()))
                {
                    //LOG_DEBUG << "\t\t[" << itr.first->ID << "] looking for var=" << depend_var;

                    auto i = std::find(itr.first->provides()->begin(), itr.first->provides()->end(), optional_var);
                    if (i != itr.first->provides()->end()) //itr provides the variable we are looking for
                    {
                        LOG_DEBUG << "\t\tAdding optional edge between " << module.first->ID << "[" <<
                                  module.first->IDnum << "] -> " << itr.first->ID << "[" << itr.first->IDnum <<
                                  "] for var=" << *i << std::endl;

                        //add the dependency from module -> itr, such that itr will come before module
                        edge e;
                        e.variable = *i;

                        //check if we need to ignore this edge as a result of user specific overrdige
                        //this will avoid a cycle is this exists.
                        bool ignore = false;
                        for (auto o : _overrides)
                        {
                            if (o.first == module.first->ID &&
                                o.second == itr.first->ID)
                            {
                                ignore = true;
                                LOG_WARNING << "Skipped adding edge between " << o.first << " and " << o.second <<
                                            " for var=" << *i << " because of user override" << std::endl;
                            }
                        }

                        if (!ignore)
                            boost::add_edge(itr.first->IDnum, module.first->IDnum, e, g);

                        //output_graph << itr.first->IDnum << "->" << module.first->IDnum << " [label=\"" << *i << "\"];" << std::endl;
                        //curr_mod_depends[*i]++; //ref count our variable

                        graphviz_vars.insert(*i);

                        module.first->set_optional_found(*i);
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
                ss << "Missing inter-module dependencies for module [" << module.first->ID << "]: " << itr.first << "\n";
                missing_depends = true;
            }

        }
        if (missing_depends)
        {
            LOG_ERROR << ss.str();
            BOOST_THROW_EXCEPTION(module_error() << errstr_info(ss.str()));
        }

        //check if our module has any met file dependenices
        if (module.first->depends_from_met()->size() == 0)
        {
            LOG_DEBUG << "Module [" << module.first->ID << "], no met file dependenices";
        } else
        {
            LOG_DEBUG << "Module [" << module.first->ID << "] has met file dependencies";
        }




        //check this modules met dependencies, bail if we are missing any.
        for (auto &depend_met_var : *(module.first->depends_from_met()))
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
        std::system("gvpr -c -f filter.gvpr -o modules.dot modules.dot.tmp");
        std::system("dot -Tpdf modules.dot -o modules.pdf");
        std::remove("modules.dot.tmp");
        std::remove("filter.gvpr");
        std::remove("modules.dot");

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
    std::system("gvpr -c -f filter.gvpr -o modules.dot modules.dot.tmp");
    std::system("dot -Tpdf modules.dot -o modules.pdf");
    std::remove("modules.dot.tmp");
    std::remove("filter.gvpr");


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

//    _ui.write_modules(s.substr(0, s.length() - 2));
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


    LOG_DEBUG << "Starting model run";


    c.tic();

    double meantime = 0;
    size_t current_ts = 0;
    _global->timestep_counter = 0; //use this to pass the timestep info to the modules for easier debugging specific timesteps
    size_t max_ts = _global->_stations.at(0)->date_timeseries().size();
    bool done = false;


        while (!done)
        {
            //ensure all the stations are at the same timestep
            boost::posix_time::ptime t;
            t = _global->_stations.at(0)->now().get_posix(); //get first stations time
            for (size_t i = 1; //on purpose to skip first station
                 i < _global->_stations.size();
                 i++)
            {
                if (t != _global->_stations.at(i)->now().get_posix())
                {
                    std::stringstream expected;
                    expected << _global->_stations.at(0)->now().get_posix();
                    std::stringstream found;
                    found << _global->_stations.at(i)->now().get_posix();
                    BOOST_THROW_EXCEPTION(forcing_timestep_mismatch()
                                          <<
                                          errstr_info("Timestep mismatch at station: " + _global->_stations.at(i)->ID()
                                                      + "\nExpected: " + expected.str()
                                                      + "\nFound: " + found.str()
                                          ));
                }
            }

            _global->_current_date = _global->_stations.at(0)->now().get_posix();


            if (!_enable_ui)
            {
                LOG_DEBUG << "Timestep: " << _global->posix_time() << "\tstep#"<<current_ts;
            }

            if(_use_netcdf)
            {
                c.tic();
                // don't use the stations variable map as it'll contain anything inserted by a filter which won't exist in the nc file
                for (auto &itr: nc.get_variable_names() )
                {
                    auto data = nc.get_var(itr, t);

                    #pragma omp parallel for
                    for (size_t y = 0; y < nc.get_ysize(); y++)
                    {
                        for (size_t x = 0; x < nc.get_xsize(); x++)
                        {
                            size_t index = x + y * nc.get_xsize();
                            auto s = _global->_stations.at(index);

                            //sanity check that we are getting the right station for this xy pair
                            if (s->ID() != std::to_string(index))
                            {
                                BOOST_THROW_EXCEPTION(
                                        forcing_error() << errstr_info("Station=" + s->ID() + ": wrong ID"));
                            }


                            double d = data[y][x];
                            s->now().set(itr, d);

                        }
                    }
                }

                LOG_DEBUG << "Done loading forcing [" << c.toc<s>() << "s]";

                c.tic();
                //do 1 step of the filters. Filters do not have depends!!
                //asume every filter is run everywhere with the same parameters
                #pragma omp parallel for
                for(size_t i = 0; i < _global->number_of_stations();i++)
                {
                    auto s = _global->stations().at(i);
                    for (auto &f : _netcdf_filters)
                    {
                        f.second->process(s);
                    }
                }

                LOG_DEBUG << "Done filters [ " << c.toc<s>() << "s]";

            }
            else
            {
                //do 1 step of the filters. Filters do not have depends!!
                #pragma omp parallel for
                for(size_t i = 0; i < _global->number_of_stations();i++)
                {
                    auto s = _global->stations().at(i);

                    //get the list of filters to run for this station
                    auto filters = _txtmet_filters[s->ID()];

                    for(auto filt : filters)
                    {
                        filt->process(s);
                    }

                }
            }

            std::stringstream ss;
            ss << _global->posix_time();
            _ui.write_timestep(ss.str());
            _ui.write_progress(int((double) current_ts / (double) max_ts * 100.0));

            c.tic();
            size_t chunks = 0;
            try
            {
                for (auto &itr : _chunked_modules)
                {
                    LOG_VERBOSE << "Working on chunk[" << chunks << "]:parallel=" <<
                                (itr.at(0)->parallel_type() == module_base::parallel::data ? "data" : "domain");

                    if (itr.at(0)->parallel_type() == module_base::parallel::data)
                    {
                        ompException oe;
                        #pragma omp parallel for
                        for (size_t i = 0; i < _mesh->size_faces(); i++)
                        {
                            auto face = _mesh->face(i);
                            if (point_mode.enable && face->_debug_name != _outputs[0].name)
                                continue;
                            oe.Run([&]
                                   {
                                       //module calls
                                       for (auto &jtr : itr)
                                       {
                                           jtr->run(face);
                                       }

                                   });
                        }
                        oe.Rethrow();

                    } else
                    {
                        //module calls for domain parallel
                        for (auto &jtr : itr)
                        {
                            ompException oe;
                            oe.Run([&]
                                   {
                                       jtr->run(_mesh);
                                   });
                            oe.Rethrow();
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
                    if(current_ts % itr.frequency == 0)
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

                                        if (jtr == output_info::mesh_outputs::vtu  )
                                        {
                                            pt::ptree &dataset = pvd.add("VTKFile.Collection.DataSet", "");
                                            dataset.add("<xmlattr>.timestep", _global->posix_time_int());
                                            dataset.add("<xmlattr>.group", "");
                                            dataset.add("<xmlattr>.part", 0);

                                            //because a full path can be provided for the base_name, we need to strip this off
                                            //to make it a relative path in the xml file.
                                            boost::filesystem::path p(base_name);

                                            dataset.add("<xmlattr>.file", p.filename().string() + ".vtu");
                                            _mesh->write_vtu(base_name + ".vtu");
                                        }

                                        if (jtr == output_info::mesh_outputs::vtp)
                                        {
                                            _mesh->write_vtp(base_name + ".vtp");
                                        }

                                    }
                                }
                            }
                        }
                    }
                }

            }


            //get data from the face into the output timeseries
            for (auto &itr : _outputs)
            {
                //only update the full timeseries
                if (itr.type == output_info::output_type::time_series)
                {
                    for (auto v : _provided_var_module)
                    {
                        auto data = itr.face->face_data(v);
                        itr.ts.at(v, current_ts) = data;
                    }
                }
            }

            if (_per_triangle_timeseries)
            {
                #pragma omp parallel for
                for (size_t i = 0; i < _mesh->size_faces(); i++)//update all the internal iterators
                {
                    if(point_mode.enable && _mesh->face(i)->_debug_name != _outputs[0].name )
                        continue;

                    _mesh->face(i)->next();
                }

            }

            //update all the stations internal iterators to point to the next time step
            for (auto &itr : _global->_stations)
            {
                if (!itr->next()) //this met station has no more met data, so doesn't matter what, we need to end now.
                {
                    done = true;
                    break;
                }
            }


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
            _ui.write_meantime(s);



            //we need it in seconds now
            if (ms)
            {
                mt /= 1000.0;
            }


            boost::posix_time::ptime pt(boost::posix_time::second_clock::local_time());
            pt = pt + boost::posix_time::seconds(mt * (max_ts - current_ts));
            _ui.write_time_estimate(boost::posix_time::to_simple_string(pt));

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


#if (BOOST_VERSION / 100 % 1000) < 56
            pt::write_xml(base_name + ".pvd",
                          pvd, std::locale(), pt::xml_writer_make_settings<char>(' ', 4));
#else
            pt::write_xml(itr.fname + ".pvd",
                          pvd, std::locale(), pt::xml_writer_settings<std::string>(' ', 4));
            break;
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

    if (_per_triangle_timeseries)
    {
#pragma omp parallel for
        for (size_t i = 0; i < _mesh->size_faces(); i++)//update all the internal iterators
        {
            auto db = _mesh->face(i)->cell_id;

            auto f = (o_path / "points") / (std::to_string(db) + "_timeseries.txt");

            _mesh->face(i)->to_file(f.string());
        }

    }

    if(_notification_script != "")
    {
        LOG_DEBUG << "Calling notification script";
        std::system(_notification_script.c_str());
    }



}

void core::end()
{
    LOG_DEBUG << "Cleaning up";
    _ui.end(); //make sure we clean up
}
