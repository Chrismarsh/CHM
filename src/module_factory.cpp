
#include "module_factory.hpp"



ModuleBase* ModuleFactory::get(std::string ID)
{
        sev_logger& lg = logger::get();
        BOOST_LOG_SEV(lg,log_level::debug) << "Requested module of type " << ID;
        
        ModuleBase* mod = NULL;
        if( ID == "test")
                mod =  new TestModule(ID);
        else
                mod =  NULL;
        
        if(mod == NULL)
            BOOST_LOG_SEV(lg,log_level::debug) << "Unable to find module";
        else
            BOOST_LOG_SEV(lg,log_level::debug) << "Succesfully instantiated module " << ID;
        
        return mod;
}
