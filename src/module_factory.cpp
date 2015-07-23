
#include "module_factory.hpp"



module_base* module_factory::get(std::string ID)
{
        LOG_DEBUG << "Module ID=" << ID;
        
        module_base* mod = NULL;
        if (ID == "solar")
            mod = new Solar(ID);
        else if (ID == "Marsh_shading_iswr")
            mod = new Marsh_shading_iswr(ID);
        else if (ID == "const_llra_ta")
            mod = new const_llra_ta(ID);
        else if (ID == "Liston_monthly_llra_ta")
            mod = new Liston_monthly_llra_ta(ID);
        else if (ID == "Liston_monthly_llra_rh")
            mod = new Liston_monthly_llra_rh(ID);
        else if (ID == "Sicart_ilwr")
            mod = new Sicart_ilwr(ID);
        else if (ID == "wind")
            mod = new wind(ID);
        else if (ID == "PenmanMonteith_evaporation")
            mod = new PenmanMonteith_evaporation(ID);
        else if (ID == "precip")
            mod = new precip(ID);
        else if (ID == "Walcek_atm_trans")
            mod = new Walcek_atm_trans(ID);


        
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
