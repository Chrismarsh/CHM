#include "scale_wind_speed.h"

scale_wind_speed::scale_wind_speed()
{

}

scale_wind_speed::~scale_wind_speed()
{



}

void scale_wind_speed::init(boost::shared_ptr<station>& station)
{
    //look at the config data to determine what we are modifying
    var = cfg.get<std::string>("variable");
    Z_F          = cfg.get<double>("Z_F"); // Get measurement height [m]
    Z_R          = Atmosphere::Z_U_R; // Reference wind speed height [m]

    // Initialize new wind speed at ref height variable
    station->add_variable("U_R");
}
void scale_wind_speed::process(boost::shared_ptr<station>& station)
{
    double U_F = station->now().get(var); // Here wind u [m/s] at Z_U
    double U_R = -9999;
    if(!is_nan(U_F))
    {
        U_R = Atmosphere::log_scale_wind(U_F, Z_F, Z_R, 0); // Assume 0 snow depth
    }

    station->now().set("U_R",U_R);

}
