
#include "module_factory.hpp"



ModuleBase* module_factory::get(std::string ID)
{
        sev_logger& lg = logger::get();
        LOG_DEBUG << "Module ID=" << ID;
        
        ModuleBase* mod = NULL;
        if( ID == "TestModule")
                mod =  new TestModule(ID);
        else
                mod =  NULL;
        
        if(mod == NULL)
        {
            BOOST_THROW_EXCEPTION(module_not_found() 
                                  << errstr_info( std::string("Module not found ") + ID)
                            );
        }
        else
        {
            return mod;
        }
}
