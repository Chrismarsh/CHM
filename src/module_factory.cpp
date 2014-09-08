
#include "module_factory.hpp"



module_base* module_factory::get(std::string ID)
{
        LOG_DEBUG << "Module ID=" << ID;
        
        module_base* mod = NULL;
        if (ID == "solar")
            mod = new Solar(ID);
        else if (ID == "shadows")
            mod = new terrain_shadow(ID);
        else if (ID == "tair_llra_const")
            mod = new tair_llra_const(ID);
        else if (ID == "tair_llra_lookup")
            mod = new tair_llra_lookup(ID);
        else if (ID == "rh_llra_var")
            mod = new rh_llra_var(ID);


        
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
