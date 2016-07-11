#include "scale_wind_speed.h"


scale_wind_speed::scale_wind_speed()
{

}
scale_wind_speed::~scale_wind_speed()
{

}
void scale_wind_speed::process(boost::shared_ptr<station> station)
{

    //look at the config data to determine what we are modifying
    std::string var = cfg.get<std::string>("variable");


    do{
        double u = station->now().get(var); // Here wind u m/s at Z_U
        double Z_U = 2.0; // cfg.get<std::double>("Z_U"); // // get input height [m] TODO: get from cfg
        double Z_U_R = 50.0; // cfg.get<std::double>("Z_U_R"); // Reference wind speed height [m] TODO: get from cfg
        
        const double Z0_SNOW = 0.01; // Snow roughness (m)
        double Ftcr;
        // Calc factor
        
        // Logrithmic wind profile assumption
        // Assumping no canopy betwen Z_U and Z_U_R

        if (Z_U > Z_U_R) { // If reference height is higher than measured height
            Ftcr = log((Z_U_R + Z0_SNOW) / Z0_SNOW) / log(Z_U / Z0_SNOW);
        } else { // If measured height is higher than reference height
            // Flip inputs, and then take 1 over Fct to get correct scalling both ways
            Ftcr = log((Z_U + Z0_SNOW) / Z0_SNOW) / log(Z_U_R / Z0_SNOW);
            Ftcr = 1 / Ftcr;
        }

        // Calc new Wind speed
        u = u * Ftcr;

        station->now().set(var,u);
    }while(station->next());
}