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
        double U_F          = station->now().get(var); // Here wind u [m/s] at Z_U
        double Z_F          = cfg.get<double>("Z_F"); // Get measurement height [m]
        double Z_R          = Atmosphere::Z_U_R; // Reference wind speed height [m]
        double U_R          = Atmosphere::log_scale_wind(U_F, Z_F, Z_R);

        station->now().set("U_R",U_R);
    }while(station->next());
}
