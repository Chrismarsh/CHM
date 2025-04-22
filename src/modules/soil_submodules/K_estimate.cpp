#include "K_estimate.hpp"

K_estimate::K_estimate(two_layer_DTO& _DTO) : DTO(_DTO)
{
    Vels = std::make_unique<Darcy_Vels>(DTO);
};

void K_estimate::run(void)
{
    Vels->init_vels();
        
    check_soil_zeros();

    Vels->calculate(); 

    set_K_values(*Vels);
           

};

void K_estimate::set_K_values(I_Darcy_Vels& Vels)
{ 
    double unit_changer = DTO.get_dt(DTO) * 1000.0; // m/s * s/step * mm/m = mm/step | m/s -> units of Vels.lateral_rechr (and others)
    double unit_changer_lateral = unit_changer / 1000.0; // dW/A from Fang et al, (2013), only applies to lateral flow
    DTO.K_rechr_to_ssr = Vels.lateral_rechr * DTO.soil_rechr_max * unit_changer_lateral;
    DTO.K_lower_to_ssr = Vels.lateral_lower * (DTO.soil_storage_max - DTO.soil_rechr_max) * unit_changer_lateral;
    DTO.K_depression_to_ssr = Vels.lateral_lower * DTO.soil_storage_max * unit_changer_lateral;
    DTO.K_depression_to_gw = Vels.vertical_depression * unit_changer;

    DTO.K_soil_to_gw = Vels.vertical_lower * unit_changer;

    DTO.K_ground_water_out = Vels.lateral_ground_water * DTO.ground_water_storage * unit_changer_lateral; // CRHM uses ground_water_storage, rather than ground_water_max in K_estimate. This was a suggestioning by a referee because ground water is saturated flow and so the "max" doesn't really make sense in this context. Remember that this model is conceptual and not a "real" soil model.
    DTO.K_detention_to_runoff = Vels.lateral_detention * DTO.detention_max * unit_changer_lateral;
};

void K_estimate::check_soil_zeros()
{
    // floating point errors can lead to some things being nearly zero but not quite
    // this function simple ensures small values are set to zero.
  
    if ( DTO.soil_rechr_storage <= 0.0000001 )
        DTO.soil_rechr_storage = 0.0;

    if ( DTO.soil_storage <= 0.0000001 )
        DTO.soil_storage = 0.0;
   
    if ( DTO.ground_water_storage <= 0.0000001 )
        DTO.ground_water_storage = 0.0;

    if ( DTO.soil_rechr_storage > DTO.soil_storage )
        DTO.soil_rechr_storage = DTO.soil_storage;

};

void I_Darcy_Vels::init_vels(void)
{
    lateral_rechr = 0.0;

    lateral_lower = 0.0;

    vertical_depression = 0.0;

    vertical_lower = 0.0;

    lateral_ground_water = 0.0;

    lateral_detention = 0.0;

};

void Darcy_Vels::calculate()
{

    if( DTO.soil_storage_max > 0.0 )
    {
        if (DTO.swe > 0.0)
        {
            set_snow();
        }
        else
            set_clear();

    }

};

void Darcy_Vels::set_snow()
{
    lateral_rechr = 0.0;
    
    lateral_lower = get_lateral_lower();
    
    vertical_depression = 0.0;
    
    vertical_lower = get_reused();
    
    lateral_ground_water = get_lateral_ground_water();
    
    if (DTO.detention_max > 0.0)
    { 
        if (DTO.snow_density > 100) // when snowcover, use Shimizu (1970) to estimate sat. hydraulic conductivity of snow
        {
            lateral_detention = get_detention_snow();
        }
        else
        {
            lateral_detention = get_detention_organic();
        }
    }
};

void Darcy_Vels::set_clear()
{
    lateral_rechr = DTO.Ksaturated_rechr * std::pow( DTO.soil_rechr_storage/DTO.soil_rechr_max, exponent) * std::tan(DTO.local_slope);
    lateral_lower = get_lateral_lower();
    
    vertical_depression = get_reused();
    
    vertical_lower = get_reused();
    
    lateral_ground_water = get_lateral_ground_water();
   
    if (DTO.detention_max > 0.0)
    { 
        lateral_detention = get_detention_organic();
    }
};

double Darcy_Vels::get_lateral_lower()
{
    return DTO.Ksaturated_lower * std::pow( (DTO.soil_storage - DTO.soil_rechr_storage) / (DTO.soil_storage_max - DTO.soil_rechr_max), exponent)  *std::tan(DTO.local_slope);
};

double Darcy_Vels::get_reused()
{
    // This Darcy Vel is used many times, so its called repeated.
    return DTO.Ksaturated_lower * std::pow( DTO.soil_storage / DTO.soil_storage_max, exponent);
};

double Darcy_Vels::get_lateral_ground_water(void)
{
    return DTO.Ksaturated_ground_water * std::tan(DTO.local_slope);
};

double Darcy_Vels::get_detention_snow(void)
{
    double Ksaturated_snow = (0.077*std::pow((DTO.snow_grain_diameter/1000),2.0)*std::exp(-7.8*(DTO.snow_density/1000)))*factor;


    return Ksaturated_snow * std::pow(DTO.detention_storage/DTO.detention_max,DTO.soil_index) * std::sin(DTO.local_slope);
}

double Darcy_Vels::get_detention_organic(void)
{
    return DTO.Ksaturated_organic * std::pow(DTO.detention_storage/DTO.detention_max,exponent_organic) * std::tan(DTO.local_slope);
};
