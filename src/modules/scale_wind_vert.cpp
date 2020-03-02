//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//


#include "scale_wind_vert.hpp"
REGISTER_MODULE_CPP(scale_wind_vert);

scale_wind_vert::scale_wind_vert(config_file cfg)
        : module_base("scale_wind_vert", parallel::data, cfg)
// Assume we are running in point mode, this gets us past the check in core. However we can later enable domain mode if we need to. We can't
// do that in this ctor as global isn't defined yet and we don't

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

void scale_wind_vert::point_scale(mesh_elem &face)
{
//    if (face->cell_local_id == 1248)
//    {
//        LOG_DEBUG << "Face found";
//    }

    // Get meteorological data for current face
    double U_R = (*face)["U_R"_s]; // Wind speed at reference height Z_R (m/s)

    // Get height info
    double Z_R = Atmosphere::Z_U_R; // Reference wind speed height [m] = 50.0 m

    double Z_CanTop = 0;

    if (!ignore_canopy && face->has_vegetation())
    {

        Z_CanTop = face->veg_attribute("CanopyHeight");
    }
    double Z_CanBot = Z_CanTop /
                      2.0; //global_param->parameters.get<double>("landcover." + std::to_string(LC) + ".TrunkHeight"); // TODO: HARDCODED until we get from obs
    //double Z_CanMid       = (Z_CanTop+Z_CanBot)/2.0; // Mid height of canopy
    double snowdepthavg = 0;
    if (has_optional("snowdepthavg"))
        snowdepthavg = (*face)["snowdepthavg"_s];

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
        (*face)["U_2m_above_srf"_s]= U_R;
        return;
    }


    /////
    // Call wind speed scaling functions to specific heights
    /////

    // If a Canopy exists
    if ( !ignore_canopy && (Z_CanTop > 0.0) && (Z_2m_above_srf <  Z_CanTop))
    {
        // Get Canopy/Surface info

        //assume we have LAI, otherwise it will cleanly bail if we don't
        double LAI = face->veg_attribute("LAI");
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

            // snowdepth is below bottom of canopy, scale bottom of canopy wind to surface
            if (Z_2m_above_srf < Z_CanBot)
                U_2m_above_srf = Atmosphere::log_scale_wind(U_CanBot, Z_CanBot, Z_2m_above_srf,
                                                            snowdepthavg);
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
    (*face)["U_2m_above_srf"_s]= U_2m_above_srf;

};

void scale_wind_vert::init(mesh& domain)
{
    //since we can optionally run this module in point or domain mode, we need to turn back on the domain parallel flag at this point
    if(!global_param->is_point_mode())
        _parallel_type =  parallel::domain;


#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
      auto face = domain->face(i);
      // Exception throwing from OpenMP needs to be here

	       auto data = face->make_module_data<d>(ID);
	       data->interp.init(interp_alg::tpspline,3,{{"reuse_LU","true"}});

    }


    ignore_canopy = cfg.get("ignore_canopy",false);

}

void scale_wind_vert::run(mesh_elem &face)
{
    point_scale(face);
}

void scale_wind_vert::run(mesh& domain)
{

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
      auto face = domain->face(i);
      // Exception throwing from OpenMP needs to be here

	       point_scale(face);

    }


#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
      // Exception throwing from OpenMP needs to be here

	       std::vector<boost::tuple<double, double, double> > u;
	       for (size_t j = 0; j < 3; j++)
	       {
		 auto neigh = face->neighbor(j);

		 if (neigh != nullptr && !neigh->_is_ghost)
		   u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(), (*neigh)["U_2m_above_srf"_s]));
	       }

	       auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());

	       if(u.size()>0)
	       {
		 double new_u =  face->get_module_data<d>(ID)->interp(u, query);
		 face->get_module_data<d>(ID)->temp_u = std::max(0.1,new_u);
	       }
	       else
	       {
		 face->get_module_data<d>(ID)->temp_u = std::max(0.1,(*face)["U_2m_above_srf"_s]);
	       }

    }


#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
      auto face = domain->face(i);
      // Exception throwing from OpenMP needs to be here

	       (*face)["U_2m_above_srf"_s]=face->get_module_data<d>(ID)->temp_u ;

    }


}
