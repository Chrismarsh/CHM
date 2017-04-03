#include "scale_wind_vert.hpp"

scale_wind_vert::scale_wind_vert(config_file cfg)
        :module_base(parallel::data)

{
    depends("U_R");

    optional("snowdepthavg");

    provides("U_2m_above_srf");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

scale_wind_vert::~scale_wind_vert()
{


}

void scale_wind_vert::run(mesh_elem &face) {

    // Bool to adjust forest edges
    bool forest_edge = cfg.get("adjust_forest_edge_wind",true);

    /////
    // Call wind speed scaling functions to specific heights
    /////

    double U_2m_above_srf = scale_wind_vert::get_U_2m_above_srf(face);

    // Adjust if cell is a forest edge
    // TODO: check if landcover exists
    if ((forest_edge) && (face->has_parameter("landcover"))) {

        // Get current cell's landcover
//        int LC = face->get_parameter("landcover");
        // 0=forest
        // 1=clearing

//        std::vector<int> n_LC = {-1, -1, -1}; // neighbor's landcover types (
        // Check if current cell's neighbors have a different landcover type
//        bool edge_f = false; // bool if current cell is on edge of any forest/clearing
        // For each neighbor
        int non_null = 1;
        for (int i = 0; i < 3; ++i) {
            auto n = face->neighbor(i); // Pointer to neighbor face
            // Check if not-null (null indicates edge cell)
            if (n != nullptr) {
                U_2m_above_srf += scale_wind_vert::get_U_2m_above_srf(n);
                non_null += 1;
//                auto n_data = n->get_module_data<snow_slide::data>(ID); // pointer to face's data
                // Calc weighting based on height diff
                // (std::max insures that if one neighbor is higher, its weight will be zero)
//                n_LC[i] = n->get_parameter("landcover");
//                if (LC != n_LC[i]) {
//                    edge_f = true;
//                }
            } // null check
        } // each neighbor
        // Take Average of neighbor cells and current cell
        U_2m_above_srf = U_2m_above_srf / non_null;
    }

    // Check that U_2m_above_srf is not too small for turbulent parameterizations (should move check there)
    if (U_2m_above_srf<0.1)
        U_2m_above_srf=0.1;

    // Save computed wind speeds
    face->set_face_data("U_2m_above_srf",U_2m_above_srf);


}

double scale_wind_vert::get_U_2m_above_srf(mesh_elem &face)
{

    // Get meteorological data for current face
    double U_R           = face->face_data("U_R"); // Wind speed at reference height Z_R (m/s)

    // Get height info
    double Z_R            = Atmosphere::Z_U_R; // Reference wind speed height [m] = 50.0 m

    double Z_CanTop       = 0;
    int LC;
    if(face->has_parameter("landcover") ) {
        LC                = face->get_parameter("landcover");
        Z_CanTop = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".CanopyHeight");
    }

    double Z_CanBot       = Z_CanTop/2.0; //global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".TrunkHeight"); // TODO: HARDCODED until we get from obs
    //double Z_CanMid       = (Z_CanTop+Z_CanBot)/2.0; // Mid height of canopy
    double snowdepthavg   = 0;
    if(has_optional("snowdepthavg"))
        snowdepthavg = face->face_data("snowdepthavg");

    // Snow depth check
    if (std::isnan(snowdepthavg)) // If it is not defined
        snowdepthavg = 0.0;

    double Z_2m_above_srf = snowdepthavg + 2.0; // (m)

    // Initialize stuff
    double U_2m_above_srf; // Wind speed 2 meters above the surface (ground or snow)
    //int LC;

    // If a Canopy exists
    if (Z_CanTop > 0.0) {
        // Get Canopy/Surface info

        //LC = face->get_parameter("landcover");
        //assume we have LAI, otherwise it will cleanly bail if we don't
        double LAI = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".LAI");
        const double alpha = LAI; // attenuation coefficient introduced by Inoue (1963) and increases with canopy density


        // Scale Z_R to Z_CanTop
        double U_CanTop = Atmosphere::log_scale_wind(U_R, Z_R, Z_CanTop, snowdepthavg);

        // Scale Z_CanTop to Z_CanBot
        double U_CanBot = Atmosphere::exp_scale_wind(U_CanTop, Z_CanTop, Z_CanBot, alpha);

        // Scale Z_CanTop to Z_CanMid
        //double U_CanMid = Atmosphere::exp_scale_wind(U_CanTop, Z_CanTop, Z_CanMid, alpha);

        // Scale Z_CanBot to Z_2m_above_srf
        if (Z_2m_above_srf < Z_CanBot) // snow depth +2 below canopy bottom
            U_2m_above_srf = Atmosphere::log_scale_wind(U_CanBot, Z_CanBot, Z_2m_above_srf,
                                                        snowdepthavg); // (U_start,Height_start,Height_end)
        else // snow depth +2 above canopy bottom
            U_2m_above_srf = Atmosphere::exp_scale_wind(U_CanTop, Z_CanTop, Z_2m_above_srf, alpha);

        // Save computed wind speeds (in case canopy exists)
        //face->set_face_data("U_CanTop",U_CanTop);
        //face->set_face_data("U_CanMid",U_CanMid);
        // No Canopy exists
    } else {
        U_2m_above_srf = Atmosphere::log_scale_wind(U_R, Z_R, Z_2m_above_srf,
                                                    snowdepthavg); // (U_start,Height_start,Height_end)
        //face->set_face_data("U_CanTop",NAN);
        //face->set_face_data("U_CanMid",NAN);
    }

    return U_2m_above_srf;
}

//void scale_wind_vert::init(mesh domain, boost::shared_ptr <global> global_param) {
//    TODO: make flag if cell is a forest edge cell (only check once here)
//}
