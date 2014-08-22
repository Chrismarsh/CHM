#include "core.h"
#include "module_base.hpp"

core::core()
{
    BOOST_LOG_FUNCTION();

    //default logging level
    _log_level = debug;

    _log_sink = boost::make_shared< text_sink >();
    {
        text_sink::locked_backend_ptr pBackend = _log_sink->locked_backend();

        boost::shared_ptr< std::ostream > pStream(&std::clog, logging::empty_deleter());
        pBackend->add_stream(pStream);


        boost::shared_ptr< std::ofstream > pStream2(new std::ofstream("CHM.log"));

        if (!pStream2->is_open())
        {
            BOOST_THROW_EXCEPTION(file_write_error()
                    << boost::errinfo_errno(errno)
                    << boost::errinfo_file_name("CHM.log")
                    );
        }

        pBackend->add_stream(pStream2);
    }

    _log_sink->set_formatter
            (
            expr::format("%1% %2% [%3%]: %4%")
            % expr::attr< boost::posix_time::ptime >("TimeStamp")
            % expr::format_named_scope("Scope",
            keywords::format = "%n:%l",
            keywords::iteration = expr::reverse,
            keywords::depth = 1)
            % expr::attr< log_level>("Severity")
            % expr::smessage
            );

    logging::core::get()->add_global_attribute("TimeStamp", attrs::local_clock());
    logging::core::get()->add_global_attribute("Scope", attrs::named_scope());

    logging::core::get()->add_sink(_log_sink);

    LOG_DEBUG << "Logger initialized. Writing to cout and CHM.log";

    #ifdef NOMATLAB
    _engine = boost::make_shared<maw::matlab_engine>();


    _engine->start();
    _engine->set_working_dir();

    LOG_DEBUG << "Matlab engine started";
#endif

    _global = boost::make_shared<global>();

    //don't just abort and die
    gsl_set_error_handler_off();
}

core::~core()
{
    LOG_DEBUG << "Terminating";
}

void core::config_debug(const json_spirit::Value& value)
{
    for (auto& jtr : value.get_obj())
    {
        const json_spirit::Pair& pair = jtr;

        if (pair.name_ == "debug_level")
        {
            std::string s = pair.value_.get_str();

            if (s == "debug")
                _log_level = debug;
            else if (s == "warning")
                _log_level = warning;
            else if (s == "error")
                _log_level = error;
            else if (s == "verbose")
                _log_level = verbose;
            else
            {
                _log_level = debug; //default to debug 
            }

            LOG_DEBUG << "Setting log severity to " << _log_level;

            _log_sink->set_filter(
                    severity >= _log_level
                    );
        }

    }
}

void core::config_modules(const json_spirit::Value& value)
{
    int modnum = 0;
    //loop over the list of requested modules
    // these are in the format "type":"ID"
    for (auto& itr : value.get_array())
    {

        const json_spirit::Value& v = itr;
        LOG_DEBUG << "Module type=" << v.get_str();

        boost::shared_ptr<module_base> m(_mfactory.get(v.get_str()));

        //internal tracking of module initialization order

        m->IDnum = modnum;
        modnum++;
        _modules.push_back(std::make_pair(m, 1)); //default to 1 for make ordering, we will set it later in determine_module_dep
    }
}

void core::config_forcing(const json_spirit::Value& value)
{
    LOG_DEBUG << "Found forcing section";
    //loop over the list of forcing data
    for (auto& jtr : value.get_obj())
    {
        const json_spirit::Pair& pair = jtr;
        const std::string& name = pair.name_;
        const json_spirit::Value& value = pair.value_;


        if (name == "station")
        {
            boost::shared_ptr<station> s = boost::make_shared<station>();
            for (auto& ktr : value.get_obj())
            {
                const json_spirit::Pair& pair = ktr;
                const std::string& name = pair.name_;


                if (name == "ID")
                {
                    s->set_ID(pair.value_.get_str());
                } else if (name == "easting")
                {
                    s->set_x(pair.value_.get_real());
                } else if (name == "northing")
                {
                    s->set_y(pair.value_.get_real());
                } else if (name == "elevation")
                {
                    s->set_z(pair.value_.get_real());
                } else if (name == "file")
                {
                    s->open(pair.value_.get_str());
                } else
                {
                    LOG_INFO << "Unknown value '" << name << "'";
                }

            }

            LOG_DEBUG << "New station created " << *s;
            _stations.push_back(s);

        } else
        {
            LOG_INFO << "Unknown forcing type " << name << ", skipping";
        }

    }
}

void core::config_meshes(const json_spirit::Value& value)
{
    LOG_DEBUG << "Found meshes section";
#ifdef MATLAB
    _mesh = boost::make_shared<triangulation>(_engine);
#else
    _mesh = boost::make_shared<triangulation>();
#endif
    for (auto& jtr : value.get_obj())
    {
        const json_spirit::Pair& pair = jtr;
        //        const std::string& name = pair.name_;
        //        const json_spirit::Value& value = pair.value_;

        if (pair.name_ == "")
        {
            BOOST_THROW_EXCEPTION(missing_value_error()
                    << errstr_info("Empty mesh name"));

            return;
        }

        std::string ID = pair.name_;
        std::string file = "";

        LOG_DEBUG << "Found mesh " << ID;

        //itr over the mesh paramters
        for (auto& ktr : pair.value_.get_obj())
        {
            const json_spirit::Pair& pair = ktr;
            //            const std::string& name = pair.name_;
            //            const json_spirit::Value& value = pair.value_;

            if (pair.name_ == "file")
            {
                file = pair.value_.get_str();
            }

        }
        _mesh->from_file(file);
        if (_mesh->size() == 0)
            BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Mesh size = 0!"));

    }
}


void core::config_matlab(const json_spirit::Value& value)
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


void core::config_output(const json_spirit::Value& value)
{
    LOG_DEBUG << "Found output section";
    //loop over the list of matlab options
    for (auto& jtr : value.get_obj())
    {
        const json_spirit::Pair& pair = jtr;
        const std::string& name = pair.name_;
        const json_spirit::Value& value = pair.value_;

        if (name == "timeseries" || name == "mesh")
        {
            output_info out;
            if (name == "timeseries")
                out.type = output_info::timeseries;
            else if (name == "mesh")
                out.type = output_info::mesh;

            LOG_DEBUG << "Found " << name;
            for (auto& ktr : value.get_obj()) //loop over all the paths
            {
                const json_spirit::Pair& pair = ktr;
                const std::string& name = pair.name_;
                const json_spirit::Value& value = pair.value_;

                if (name == "northing")
                    out.northing = value.get_real();
                else if (name == "easting")
                    out.easting = value.get_real();
                else if (name == "variables")
                {
                    for (auto& itr_vars : value.get_array())
                    {
                        const json_spirit::Value& v = itr_vars;
                        LOG_VERBOSE << "Output variable: " << v.get_str();
                        out.variables.push_back(v.get_str());
                    }
                } else if (name == "file")
                    out.out_file = value.get_str();
                else if (name == "plot")
                    out.plot = value.get_bool();
            }

            if (name == "timeseries")
            {
                out.face = _mesh->locate_face(out.easting, out.northing);
                if (out.face == NULL)
                {
                    LOG_WARNING << "Requested an output point that is not in the triangulation domain, skipping";
                } else
                {
                    _outputs.push_back(out);
                }
            } else
            {
                _outputs.push_back(out);
            }

        } else
        {
            LOG_WARNING << "Unknown output type: " << name;
        }
    }
}

void core::config_global(const json_spirit::Value& value)
{
    LOG_DEBUG << "Found global section";
    for (auto& jtr : value.get_obj())
    {
        const json_spirit::Pair& pair = jtr;
        const std::string& name = pair.name_;
        const json_spirit::Value& value = pair.value_;

        if (name == "latitude")
        {
            _global->_lat = value.get_real();
        } else if (name == "longitude")
        {
            _global->_lon = value.get_real();
        } else if (name == "UTC_offset")
        {
            _global->_utc_offset = value.get_int();
        }
    }
}

void core::read_config_file(std::string file)
{
    BOOST_LOG_FUNCTION();
    std::ifstream is(file);
    if (is.fail())
    {
        BOOST_THROW_EXCEPTION(file_read_error()
                << boost::errinfo_errno(errno)
                << boost::errinfo_file_name(file)
                );

    }

    LOG_DEBUG << "Reading configuration file " << file;
    //holds the top level object { ... }
    json_spirit::Value value;

    json_spirit::read(is, value);

    //get the top-level 
    const json_spirit::Object& top_level = value.get_obj();


    //loop over the top-level options
    for (auto& itr : top_level)
    {
        const json_spirit::Pair& pair = itr;
        const std::string& name = pair.name_;
        const json_spirit::Value& value = pair.value_;

        //check for debug specification
        if (name == "debug")
        {
            config_debug(value);
        } else if (name == "modules")
        {
            config_modules(value);

        } else if (name == "forcing")
        {
            config_forcing(value);
        } else if (name == "meshes")
        {
            config_meshes(value);
        } else if (name == "matlab")
        {
            #ifdef NOMATLAB
            config_matlab(value);
#endif
        } else if (name == "output")
        {
            config_output(value);
        } else if (name == "global")
        {
            config_global(value);
        } else
        {
            const json_spirit::Pair& pair = itr;
            const std::string& name = pair.name_;
            LOG_INFO << "Unknown section '" << name << "', skipping";
        }
    }

    //    if (found_modules == false)
    //    {
    //        BOOST_THROW_EXCEPTION(no_modules_defined()
    //                << errstr_info(std::string("No module section found in ") + file)
    //                );
    //    }


    LOG_DEBUG << "Finished initialization";

    LOG_DEBUG << "Init variables mapping";
    _global->_variables.init_from_file("Un-init path");

    LOG_DEBUG << "Determining module dependencies";
    _determine_module_dep();

}

void core::_determine_module_dep()
{
    int size = _modules.size();

    Graph g(size);
    //    std::vector<Edge> edges;

    //loop through each module
    for (auto& module : _modules)
    {
        //Generate a debugging list of all variables, culling duplicates
        _module_provided_variable_list.insert(module.first->provides()->begin(), module.first->provides()->end());
        //look at all other modules
        if (module.first->depends()->size() == 0)
        {
            LOG_DEBUG << "Checking [" << module.first->ID << "], No dependenices";
        } else
        {
            LOG_DEBUG << "Checking [" << module.first->ID << "] against...";
        }

        for (auto& itr : _modules)
        {
            //don't check against our module
            if (module.first->ID.compare(itr.first->ID) != 0)
            {
                //loop through each required variable of our current module
                for (auto& depend_var : *(module.first->depends()))
                {
                    LOG_DEBUG << "Module=" << itr.first->ID << " looking for var=" << depend_var;

                    auto i = std::find(itr.first->provides()->begin(), itr.first->provides()->end(), depend_var);
                    if (i != itr.first->provides()->end()) //itr provides the variable we are looking for
                    {
                        LOG_DEBUG << "Adding edge between " << module.first->ID << "[" << module.first->IDnum << "] -> " << itr.first->ID << "[" << itr.first->IDnum << "] for var=" << *i << std::endl;

                        //add the dependency from module -> itr, such that itr will come before module
                        boost::add_edge(module.first->IDnum, itr.first->IDnum, g);
                    }
                }
            }
        }
    }

    std::list<Vertex> make_order;

    boost::topological_sort(g, std::front_inserter(make_order));

    for (auto i = make_order.rbegin(); i != make_order.rend(); ++i)
    {
        _modules.at(*i).second = *i;
    }

    //sort descending
    std::sort(_modules.begin(), _modules.end(),
            [](const std::pair<module, size_t>& a, const std::pair<module, size_t>& b)->bool
            {
                return a.second < a.second;
            });


    std::stringstream ss;
    for (auto itr : _modules)
    {
        ss << itr.first->ID << "->";

    }
    std::string s = ss.str();
    LOG_DEBUG << "Build order: " << s.substr(0, s.length() - 2);

    // Parallel compilation ordering
    std::vector<int> time(size, 0);
    for (auto i = make_order.begin(); i != make_order.end(); ++i)
    {
        // Walk through the in_edges an calculate the maximum time.
        if (boost::in_degree(*i, g) > 0)
        {
            Graph::in_edge_iterator j, j_end;
            int maxdist = 0;
            // Through the order from topological sort, we are sure that every 
            // time we are using here is already initialized.
            for (boost::tie(j, j_end) = boost::in_edges(*i, g); j != j_end; ++j)
                maxdist = (std::max)(time[boost::source(*j, g)], maxdist);
            time[*i] = maxdist + 1;
        }
    }
    boost::graph_traits<Graph>::vertex_iterator i, iend;
    for (boost::tie(i, iend) = boost::vertices(g); i != iend; ++i)
    {
        LOG_DEBUG << "time_slot[" << time[*i] << "] = " << _modules.at(*i).first->ID << std::endl;
    }

    for (size_t i = 0; i < _stations.size(); i++)
    {
        auto vars = _stations.at(i)->list_variables();
        _module_provided_variable_list.insert(vars.begin(), vars.end());
    }

    LOG_DEBUG << "List of all provided variables: ";
    std::string modlist = "";
    for (auto itr : _module_provided_variable_list)
    {
        s += itr + " ";
    }
    LOG_DEBUG << modlist;

    LOG_DEBUG << "Initializing and allocating memory for timeseries";

    if (_stations.size() == 0)
        BOOST_THROW_EXCEPTION(forcing_no_stations() << errstr_info("no stations"));

    //#pragma omp parallel for
    for (triangulation::Finite_faces_iterator fit = _mesh->finite_faces_begin(); fit != _mesh->finite_faces_end(); ++fit)
    {

        auto date = _stations.at(0)->get_date_timeseries();
        auto size = _stations.at(0)->get_timeseries_length();
        Delaunay::Face_handle face = fit;
        face->init_time_series(_module_provided_variable_list, /*list of all the variables that are provided by met files or modules*/
                date, /*take the first station, later checks ensure all the stations' timeseries match*/
                size); /*length of all the vectors to initialize*/
    }

}

void core::run()
{
    LOG_DEBUG << "Entering main loop";


    timer c;
    c.tic();

    interp_t_air interp;
    interp_rh irh;

    //do the interpolation first
    bool done = false;
    while (!done)
    {
        //ensure all the stations are at the same timestep
        boost::posix_time::ptime t;
        t = _stations.at(0)->now().get_posix(); //get first stations time
        for (size_t i = 1; //on purpose to skip first station
                i < _stations.size();
                i++)
        {
            if (t != _stations.at(i)->now().get_posix())
            {
                BOOST_THROW_EXCEPTION(forcing_timestep_mismatch()
                        << errstr_info("Timestep mismatch at station: " + _stations.at(i)->get_ID()));
            }
        }

        //start date
        _global->_current_date = _stations.at(0)->now().get_posix();

        //calculate global, e.g. solar position
        _global->update();


        LOG_DEBUG << "Interpolating at timestep: " << _global->posix_time();

        //iterate over all the mesh elements
        for (triangulation::Finite_faces_iterator fit = _mesh->finite_faces_begin(); fit != _mesh->finite_faces_end(); ++fit)
        {
            //interpolate the station data to the current element
            triangulation::Face_handle face = fit;
            interp("LLRA_var", face, _stations, _global);
            irh("LLRA_rh_var", face, _stations, _global);

            //this triangle needs to be advanced to the next timestep
            fit->next();

        }

        //update all the stations internal iterators to point to the next time step
        for (auto& itr : _stations)
        {
            if (!itr->next()) //
                done = true;
        }

    }

    //reset all the internal iterators
    for (auto& itr : _stations)
    {
        itr->reset_itrs();
    }

    //reset the iterators for all mesh timeseries
    //#pragma omp parallel for
    for (triangulation::Finite_faces_iterator fit = _mesh->finite_faces_begin(); fit != _mesh->finite_faces_end(); ++fit)
    {
        //current mesh element
        fit->reset_to_begining();
    }


    done = false;
    while (!done)
    {

        _global->_current_date = _stations.at(0)->now().get_posix();
        _global->update();


        LOG_DEBUG << "Timestep: " << _global->posix_time();

        //iterate over all the mesh elements
        //#pragma omp parallel for
        for (triangulation::Finite_faces_iterator fit = _mesh->finite_faces_begin(); fit != _mesh->finite_faces_end(); ++fit)
        {
            //module calls
            for (auto& itr : _modules)
            {
                triangulation::Face_handle face = fit;
                if (itr.first->parallel_type() == module_base::parallel::data)
                    itr.first->run(face, _global);

            }
        }

        //module calls for domain parallel
        for (auto& itr : _modules)
        {
            if (itr.first->parallel_type() == module_base::parallel::domain)
                itr.first->run(_mesh, _global);

        }

        //update all the stations internal iterators to point to the next time step
        for (auto& itr : _stations)
        {
            if (!itr->next()) //
                done = true;
        }
    }
    double elapsed = c.toc();
    LOG_DEBUG << "Took " << elapsed << "s";

#ifdef NOMATLAB
    for (auto& itr : _outputs)
    {
        if (itr.plot)
        {
            //loop through all the request variables
            for (auto& jtr : itr.variables)
            {
                if (itr.type == output_info::mesh)
                {
                    _mesh->plot(jtr);
                } else if (itr.type == output_info::timeseries)
                {
                    _mesh->plot_time_series(itr.easting, itr.northing, jtr);
                }
            }
        }
        if (itr.out_file != "")
        {
            _mesh->to_file(itr.face, itr.out_file);
        }
    }
#endif
}