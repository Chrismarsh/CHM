#include "core.h"
#include "module_base.hpp"

core::core()
{
    BOOST_LOG_FUNCTION();
    _module_config = "";

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

    _engine = boost::make_shared<maw::matlab_engine>();

    
    _engine->start();
    _engine->set_working_dir();
   
    LOG_DEBUG << "Matlab engine started";
    
    _engine->add_dir_to_path("/home/chris/Documents/PhD/code/CHM/src/matlab_support");
    
}



core::~core()
{

}

bool core::is_debug()
{
    return _is_debug;
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
    bool found_modules = false;

    //loop over the top-level options
    for (auto& itr : top_level)
    {
        const json_spirit::Pair& pair = itr;
        const std::string& name = pair.name_;
        const json_spirit::Value& value = pair.value_;


        //check for debug specification
        if (name == "debug")
        {
            for (auto& jtr : value.get_obj())
            {
                const json_spirit::Pair& pair = jtr;

                if (pair.name_ == "debug_level")
                {
                    std::string s = pair.value_.get_str();

                    if (s == "debug")
                        _log_level = log_level::debug;
                    else if (s == "warning")
                        _log_level = log_level::warning;
                    else if (s == "error")
                        _log_level = log_level::error;
                    else
                    {
                        _log_level = log_level::debug; //default to debug 
                    }

                    LOG_DEBUG << "Setting log severity to " << _log_level;

                    _log_sink->set_filter(
                            severity >= _log_level
                            );



                }

            }
        } else if (name == "modules")
        {
            const json_spirit::Pair& pair = itr;
            const std::string& name = pair.name_;
            const json_spirit::Value& value = pair.value_;
            found_modules = true; //found the module section

            //loop over the list of requested modules
            // these are in the format "type":"ID"
            for (auto& jtr : value.get_obj())
            {
                const json_spirit::Pair& pair = jtr;
                const std::string& name = pair.name_;
                const std::string& ID = pair.value_.get_str();

                LOG_DEBUG << "Module type=" << name << " ID=" << ID;
                ModuleBase* lol = _mfactory.get(ID);
                delete lol;

            }

        } else if (name == "forcing")
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
                            LOG_DEBUG << "Unknown value '" << name << "'";
                        }

                    }

                    LOG_DEBUG << "New station created " << *s;
                    _stations.push_back(s);

                } else
                {
                    LOG_DEBUG << "Unknown forcing type " << name << ", skipping";
                }


            }
        } else if (name == "meshes")
        {
            LOG_DEBUG << "Found meshes section";
            _mesh = boost::make_shared<mesh>(_engine);
            for (auto& jtr : value.get_obj())
            {
                const json_spirit::Pair& pair = jtr;
                const std::string& name = pair.name_;
                const json_spirit::Value& value = pair.value_;

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
                for (auto& ktr : value.get_obj())
                {
                    const json_spirit::Pair& pair = ktr;
                    const std::string& name = pair.name_;
                    const json_spirit::Value& value = pair.value_;

                    if (pair.name_ == "file")
                    {
                        file = pair.value_.get_str();
                    }
                    
                }
                _mesh->add_mesh(file, ID);

            }
        } else
        {
            const json_spirit::Pair& pair = itr;
            const std::string& name = pair.name_;
            LOG_DEBUG << "Unknown section '" << name << "', skipping";
        }
    }

    if (found_modules == false)
    {
        BOOST_THROW_EXCEPTION(no_modules_defined()
                << errstr_info(std::string("No module section found in ") + file)
                );
    }


    LOG_DEBUG << "Finished initialization";

}

void core::run()
{
    
    
    //iterate over all the mesh elements
//    for(size_t i=0;
//            //this assumes that all the meshes are the same size
//            i<_mesh->size(); // TODO: Add a check to ensure all meshes are the same size
//            i++)
//    {
//        
//        //interpolate the forcing data over the mesh
//        
//    }
    _mesh->plot("DEM");
    _engine->evaluate("save lol.mat");
    
    
    
    

    //run all selected algorithms

}