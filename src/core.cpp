#include "core.h"

core::core()
{
    _module_config = "";
    
    _log_level = debug;
    
    _log_sink = boost::make_shared< text_sink >();
    { 
       text_sink::locked_backend_ptr pBackend = _log_sink->locked_backend();

       boost::shared_ptr< std::ostream > pStream(&std::clog, logging::empty_deleter());
       pBackend->add_stream(pStream);


       boost::shared_ptr< std::ofstream > pStream2(new std::ofstream("CHM.log"));
       assert(pStream2->is_open());
       pBackend->add_stream(pStream2);
    }

   _log_sink->set_formatter
   (
       expr::format("[%1%] [%2%]: %3%")
           % expr::attr< boost::posix_time::ptime  >("TimeStamp")
           % expr::attr< log_level>("Severity")
           % expr::smessage
   );

    logging::core::get()->add_global_attribute("TimeStamp", attrs::local_clock());
    
    logging::core::get()->add_sink(_log_sink);
    
}

core::~core()
{
    
}

bool core::is_debug()
{
  return _is_debug;
}

void core::read_module_file(std::string file)
{
    std::ifstream is( file );
    
    if (is.fail())
    {
        BOOST_THROW_EXCEPTION(file_read_error() 
                                << boost::errinfo_errno(errno)
                                << boost::errinfo_file_name(file)
                            );
        
    }

    //holds the top level object { ... }
    json_spirit::Value value;

    json_spirit::read( is, value );
 
    //get the top-level 
    const json_spirit::Object& top_level = value.get_obj();

    
    //loop over the top-level options
    for(auto& itr : top_level )
    {
        const json_spirit::Pair& pair = itr;
        const std::string& name  = pair.name_;
        const json_spirit::Value&  value = pair.value_;
        
        std::cout << name << std::endl;
        
    }
}

void core::read_config_file( std::string file )
{
    std::ifstream is( file );
    if (is.fail())
    {
        BOOST_THROW_EXCEPTION(file_read_error() 
                                << boost::errinfo_errno(errno) 
                                << boost::errinfo_file_name(file)
                                << boost::errinfo_type_info_name("LOL")
                            );

    }


    //holds the top level object { ... }
    json_spirit::Value value;

    json_spirit::read( is, value );
 
    //get the top-level 
    const json_spirit::Object& top_level = value.get_obj();

    
    //loop over the top-level options
    for(auto& itr : top_level )
    {
        const json_spirit::Pair& pair = itr;
        const std::string& name  = pair.name_;
        const json_spirit::Value&  value = pair.value_;
	
	std::cout << name <<  std::endl;
	
	//check for debug specification
	if (name == "debug")
	{
	  for(auto& jtr : value.get_obj())
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
                  _log_level  = log_level::error;
                else 
                {
                  _log_level = log_level::debug; //default to debug 
                }
	                      
                _log_sink->set_filter(
                   severity >= _log_level
                );
                
		
	    }

	  }
	}
        else if(name == "modules")
        {
            _module_config = pair.value_.get_str();
            BOOST_LOG_SEV(logger::get(),log_level::debug) << "Found module config file at " << _module_config;
            read_module_file(_module_config);
            
        }
    }
    
    
    BOOST_LOG_SEV(logger::get(),log_level::debug) << "Successfully opened the configuration file: " << file;
    BOOST_LOG_SEV(logger::get(),log_level::debug) << "Writting logs to stdout and CHM.log at level: " << _log_level;
}
