#include "scale_wind_vert.hpp"

scale_wind_vert::scale_wind_vert(config_file cfg)
        :module_base(parallel::data)

{
    //depends("U_R"); // implicite depends on the filter scale_wind_speed

    provides("U_R");
    provides("U_CanTop");
    //provides("U_CanMid");
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
    int LC                = face->get_parameter("landcover"); // TODO: Have to call this here because Simple_Canopy::init does it (how to only do once?)
    double Z_CanTop       = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".CanopyHeight");
    double Z_CanBot       = Z_CanTop/2.0; //global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".TrunkHeight"); // TODO: HARDCODED. Get from params
    double snowdepthavg   = face->face_data("snowdepthavg");

    // Snow depth check
    if (is_nan(snowdepthavg)) // If it is not defined TODO: current hack, should be required to be initialized in mesher
        snowdepthavg = 0;

    double Z_2m_above_srf = snowdepthavg + 2.0; // (m)

    // Get Canopy/Surface info
    double LAI            = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".LAI");

    // Scale Z_R to Z_CanTop
    double U_CanTop       = Atmosphere::log_scale_wind(U_R, Z_R, Z_CanTop);

    // Scale Z_CanTop to Z_CanMid

    // Scale Z_CanTop to Z_CanBot
    double U_CanBot  = Atmosphere::log_scale_wind(U_CanTop, Z_CanTop, Z_CanBot); // (U_start,Height_start,Height_end)
    // double U_CanMid  = Atmosphere::exp_scale_wind(U_CanTop, Z_R, Z_2m_above_srf);

    // Scale Z_CanBot to Z_2m_above_srf TODO: Might need check that snowdepth is less than truck height
    double U_2m_above_srf  = Atmosphere::log_scale_wind(U_CanBot, Z_CanBot, Z_2m_above_srf); // (U_start,Height_start,Height_end)


    // Save computed wind speeds
    face->set_face_data("U_R",U_R);
    face->set_face_data("U_CanTop",U_CanTop);
    face->set_face_data("U_2m_above_srf",U_2m_above_srf);

}

void scale_wind_vert::init(mesh domain, boost::shared_ptr <global> global_param) {

}
