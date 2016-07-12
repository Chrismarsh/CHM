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
        double u              = station->now().get(var); // Here wind u [m/s] at Z_U
        double Z_U            = cfg.get<double>("Z_U"); // Get measurement height [m]
        double Z_U_R          = Atmopshere::Z_U_R; // Reference wind speed height [m] = 50.0 m

        double u_new          = log_scale_wind(u, Z_U, Z_U_R);

        station->now().set(var,u_new);
    }while(station->next());
}