
#include "scale_wind_vert.hpp"

scale_wind_vert::scale_wind_vert(config_file cfg)
        :module_base(parallel::domain)

{
    depends("U_R");

    optional("snowdepthavg");

    //provides("U_CanTop"); // Possible output, but commented out for speed
    //provides("U_CanMid");
    provides("U_2m_above_srf");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

scale_wind_vert::~scale_wind_vert()
{


}

void scale_wind_vert::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto data = face->make_module_data<d>(ID);
        data->interp.init(interp_alg::tpspline);
    }

}

void scale_wind_vert::run(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        // Get meteorological data for current face
        double U_R = face->face_data("U_R"); // Wind speed at reference height Z_R (m/s)

        // Get height info
        double Z_R = Atmosphere::Z_U_R; // Reference wind speed height [m] = 50.0 m

        double Z_CanTop = 0;

        if (face->has_parameter("landcover"))
        {
            int LC = face->get_parameter("landcover");
            Z_CanTop = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".CanopyHeight");
        }
        double Z_CanBot = Z_CanTop /
                          2.0; //global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".TrunkHeight"); // TODO: HARDCODED until we get from obs
        //double Z_CanMid       = (Z_CanTop+Z_CanBot)/2.0; // Mid height of canopy
        double snowdepthavg = 0;
        if (has_optional("snowdepthavg"))
            snowdepthavg = face->face_data("snowdepthavg");

        // Snow depth check
        if (is_nan(snowdepthavg)) // If it is not defined
            snowdepthavg = 0.0;

        double Z_2m_above_srf = snowdepthavg + 2.0; // (m)

        // Initialize stuff
        double U_2m_above_srf; // Wind speed 2 meters above the surface (ground or snow)

        // Check if snowdepth is above U_R (avalanche gone crazy case)
        // This is outside of log_scale_wind()'s range, so make assumption that wind is constant above U_R
        if (Z_2m_above_srf >= Z_R)
        {
            face->set_face_data("U_2m_above_srf", U_R);
            continue;
        }


        /////
        // Call wind speed scaling functions to specific heights
        /////

        // If a Canopy exists
        if (Z_CanTop > 0.0)
        {
            // Get Canopy/Surface info

            int LC = face->get_parameter("landcover");
            //assume we have LAI, otherwise it will cleanly bail if we don't
            double LAI = global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".LAI");
            const double alpha = LAI; // attenuation coefficient introduced by Inoue (1963) and increases with canopy density

            // If snowdepth is below the Canopy Top
            if (snowdepthavg < Z_CanTop)
            {
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
            } else
            {
                // Scale Z_R to snowdepth (which is above canopy height)
                U_2m_above_srf = Atmosphere::log_scale_wind(U_R, Z_R, Z_2m_above_srf, snowdepthavg);
            }

            // No Canopy exists
        } else
        {
            U_2m_above_srf = Atmosphere::log_scale_wind(U_R, Z_R, Z_2m_above_srf,
                                                        snowdepthavg); // (U_start,Height_start,Height_end)
        }

        // Check that U_2m_above_srf is not too small for turbulent parameterizations (should move check there)
        U_2m_above_srf = std::max(0.1, U_2m_above_srf);
        face->set_face_data("U_2m_above_srf",U_2m_above_srf);
        face->make_module_data<d>(ID);

    }
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);
        std::vector<boost::tuple<double, double, double> > u;
        for (size_t j = 0; j < 3; j++)
        {
            auto neigh = face->neighbor(j);

            if (neigh != nullptr)
                u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(), neigh->face_data("U_2m_above_srf")));


        }

        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());

        interpolation interp(interp_alg::tpspline);
        double new_u = interp(u, query);
        face->get_module_data<d>(ID)->temp_u = new_u;
    }

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->set_face_data("U_2m_above_srf",face->get_module_data<d>(ID)->temp_u );
//       face->set_face_data("U_R",x[i] );

    }

    // Save computed wind speeds
//    face->set_face_data("U_2m_above_srf",U_2m_above_srf);
}
