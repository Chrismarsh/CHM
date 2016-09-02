#include "scale_wind_vert.hpp"

scale_wind_vert::scale_wind_vert(config_file cfg)
        :module_base(parallel::data)

{
    //depends("U_R"); // implicit depends on the filter scale_wind_speed

    provides("U_R");
    provides("U_CanTop");
    provides("U_CanMid");
    provides("U_2m_above_srf");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

scale_wind_vert::~scale_wind_vert()
{


}

void scale_wind_vert::run(mesh_elem &face, boost::shared_ptr <global> global_param) {

    // TODO: Add option if the filter::scale_wind_speed was not used

    // Get meteorological data for current face
    double U_R           = face->face_data("U_R"); // Wind speed at reference height Z_R (m/s)

    // Get height info
    double Z_R            = Atmosphere::Z_U_R; // Reference wind speed height [m] = 50.0 m
    int LC                = face->get_parameter("landcover"); // TODO: Move to init
    double Z_CanTop       = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".CanopyHeight");
    double Z_CanBot       = Z_CanTop/2.0; //global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".TrunkHeight"); // TODO: HARDCODED until we get from obs
    double Z_CanMid       = (Z_CanTop+Z_CanBot)/2.0; // Mid height of canopy
    double snowdepthavg   = face->face_data("snowdepthavg");
    // Snow depth check
    if (is_nan(snowdepthavg)) // If it is not defined
        snowdepthavg = 0;
    double Z_2m_above_srf = snowdepthavg + 2.0; // (m)

    // Get Canopy/Surface info
    double LAI            = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".LAI"); // TODO: Move to init
    const double alpha    = LAI; // attenuation coefficient introduced by Inoue (1963) and increases with canopy density

    // Initialize stuff
    double U_2m_above_srf; // Wind speed 2 meters above the surface (ground or snow)

    /////
    // Call wind speed scaling functions to specific heights
    /////

    // If a Canopy exists
    if (Z_CanTop>0) {
        // Scale Z_R to Z_CanTop
        double U_CanTop = Atmosphere::log_scale_wind(U_R, Z_R, Z_CanTop, snowdepthavg);

        // Scale Z_CanTop to Z_CanBot
        double U_CanBot = Atmosphere::exp_scale_wind(U_CanTop, Z_CanTop, Z_CanBot, alpha);

        // Scale Z_CanTop to Z_CanMid
        double U_CanMid = Atmosphere::exp_scale_wind(U_CanTop, Z_CanTop, Z_CanMid, alpha);

        // Scale Z_CanBot to Z_2m_above_srf
        if (Z_2m_above_srf < Z_CanBot) // snow depth +2 below canopy bottom
            U_2m_above_srf = Atmosphere::log_scale_wind(U_CanBot, Z_CanBot, Z_2m_above_srf, snowdepthavg); // (U_start,Height_start,Height_end)
        else // snow depth +2 above canopy bottom
            U_2m_above_srf = Atmosphere::exp_scale_wind(U_CanTop, Z_CanTop, Z_2m_above_srf, alpha);

        // Save computed wind speeds (in case canopy exists)
        face->set_face_data("U_CanTop",U_CanTop);
        face->set_face_data("U_CanMid",U_CanMid);
    // No Canopy exists
    } else {
        U_2m_above_srf = Atmosphere::log_scale_wind(U_R, Z_R, Z_2m_above_srf,snowdepthavg); // (U_start,Height_start,Height_end)
    }

    // Save computed wind speeds
    face->set_face_data("U_R",U_R);
    face->set_face_data("U_2m_above_srf",U_2m_above_srf);

}

void scale_wind_vert::init(mesh domain, boost::shared_ptr <global> global_param) {

}
