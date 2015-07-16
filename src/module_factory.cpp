
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
        else if (ID == "atm_trans_annandale")
            mod = new atm_trans_annandale(ID);
        else if (ID == "longwave_sicart")
            mod = new longwave_sicart(ID);
        else if (ID == "wind")
            mod = new wind(ID);
        else if (ID == "evap_penman_monteith")
            mod = new evap_penman_monteith(ID);
        else if (ID == "precip")
            mod = new precip(ID);
        else if (ID == "leaky_bucket")
            mod = new leaky_bucket(ID);


        
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
