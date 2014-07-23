
#include "module_factory.hpp"



module_base* module_factory::get(std::string ID)
{
        sev_logger& lg = logger::get();
        LOG_DEBUG << "Module ID=" << ID;
        
        module_base* mod = NULL;
        if (ID == "Solar")
            mod = new Solar(ID);
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
