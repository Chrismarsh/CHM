
#include "module_factory.hpp"



module_base* module_factory::get(std::string ID, pt::ptree& config)
{
    LOG_DEBUG << "Module ID=" << ID;
    
    module_base* mod = nullptr;
    
    if (ID == "solar")
        mod = new slope_iswr();
    else if (ID == "Marsh_shading_iswr")
        mod = new Marsh_shading_iswr();
    else if (ID == "const_llra_ta")
        mod = new const_llra_ta();
    else if (ID == "Liston_monthly_llra_ta")
        mod = new Liston_monthly_llra_ta();
    else if (ID == "Liston_monthly_llra_rh")
        mod = new Liston_monthly_llra_rh();
    else if (ID == "Sicart_ilwr")
        mod = new Sicart_ilwr();
    else if (ID == "Liston_wind")
        mod = new Liston_wind();
    else if (ID == "PenmanMonteith_evaporation")
        mod = new PenmanMonteith_evaporation();
    else if (ID == "precip")
        mod = new precip();
    else if (ID == "Walcek_cloud")
        mod = new Walcek_cloud();
    else if (ID == "Harder_precip_phase")
        mod = new Harder_precip_phase();
    else if (ID == "Burridge_iswr")
        mod = new Burridge_iswr();
    else if (ID == "Iqbal_iswr")
        mod = new Iqbal_iswr();
    else if (ID == "iswr_from_obs")
        mod = new iswr_from_obs();
    else if (ID == "Dodson_NSA_ta")
        mod = new Dodson_NSA_ta();

    if(mod == nullptr)
    {
        BOOST_THROW_EXCEPTION(module_not_found() 
                              << errstr_info( std::string("Module not found ") + ID)
                        );
    }
    
    mod->ID = ID;
    mod->cfg = config;
    
    return mod;
       


}
