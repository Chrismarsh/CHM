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

#include "Ninja_wind.hpp"
REGISTER_MODULE_CPP(Ninja_wind);

Ninja_wind::Ninja_wind(config_file cfg)
        : module_base("Ninja_wind", parallel::domain, cfg)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("Ninja_speed");
    provides("Ninja_speed2");
    provides("vw_dir");
    provides("Ninja_u");
    provides("Ninja_v");

    provides("interp_zonal_u");
    provides("interp_zonal_v");

    provides("lookup_d");

    provides("vw_dir_orig");

    ninja_average = cfg.get("ninja_average",true);


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void Ninja_wind::init(mesh domain)
{
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
        d->interp_smoothing.init(interp_alg::tpspline,3,{ {"reuse_LU","true"}});
    }
}


void Ninja_wind::run(mesh domain)
{
    double H_forc = 40.0;  // Reference height for GEM forcing and WindNinja wind field library

    double Max_spdup = 1.5;  // Maximal value of crest speedup
    double Min_spdup = 0.1;  // Minimal value of crest speedup

        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);

            std::vector<boost::tuple<double, double, double> > u;
            std::vector<boost::tuple<double, double, double> > v;
            for (auto &s : global_param->get_stations(face->get_x(), face->get_y()))
            {
                if (is_nan(s->get("U_R")) || is_nan(s->get("vw_dir")))
                    continue;

                double theta = s->get("vw_dir") * M_PI / 180.;

                double W = s->get("U_R");
                W = std::max(W, 0.1);

                W = Atmosphere::log_scale_wind(W,
                                               Atmosphere::Z_U_R,  // UR is at our reference height
                                               H_forc,  //  Reference height for GEM forcing and WindNinja wind field library
                                               0); // no canopy, no snow, but uses a snow roughness

                double zonal_u = -W * sin(theta);
                double zonal_v = -W * cos(theta);

                u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
                v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
            }

            //http://mst.nerc.ac.uk/wind_vect_convs.html

            // get an interpolated zonal U,V at our face
            auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
            double zonal_u = face->get_module_data<data>(ID)->interp(u, query);
            double zonal_v = face->get_module_data<data>(ID)->interp(v, query);

            face->set_face_data("interp_zonal_u", zonal_u);
            face->set_face_data("interp_zonal_v", zonal_v);

            //Get back the interpolated wind direction
            // -- not sure if there is a better way to do this, but at least a first order to getting the right direction
            double theta = math::gis::zonal2dir(zonal_u, zonal_v);

            double theta_orig = theta;

            double U = 0.;
            double V = 0.;
            double W_transf = 0.;
            

            if(!ninja_average)  // No Linear interpolation between the closest 2 wind fields from the library
            {
                // Use this wind dir to figure out which wind field from the library we need
                // Wind field are available each 15 deg.
                int d = int(theta * 180.0 / M_PI / 15.);
                if (d == 0) d = 24;
                face->set_face_data("lookup_d", d);

                // get the transfert function and associated wind component for the interpolated wind direction
                 W_transf = face->get_parameter("Ninja" + std::to_string(d));   // transfert function
                 U = face->get_parameter("Ninja" + std::to_string(d) + "_U");  // zonal component
                 V = face->get_parameter("Ninja" + std::to_string(d) + "_V");  // meridional component

           }else // Linear interpolation between the closest 2 wind fields from the library
           {
 
                // Use this wind dir to figure out which wind fields from the library we need
                // Wind fields are available each 15 deg.        
                int d1 = int(theta * 180.0 / M_PI / 15.);
                double theta1 = d1 * 15.0 * M_PI / 180.0;
                if (d1 == 0) d1 = 24;

                int d2 = int((theta * 180.0 / M_PI + 15.) / 15.);
                double theta2 = d2 * 15.0 *M_PI / 180.0;
                if (d2 == 0) d2 = 24;

                double d = d1*(theta2-theta)/(theta2-theta1)+d2*(theta-theta1)/(theta2-theta1);
                face->set_face_data("lookup_d", d);

                // get the transfert function and associated wind component for the interpolated wind direction
                double W_transf1 = face->get_parameter("Ninja" + std::to_string(d1));   // transfert function
                double U_lib1 = face->get_parameter("Ninja" + std::to_string(d1) + "_U");  // zonal component
                double V_lib1 = face->get_parameter("Ninja" + std::to_string(d1) + "_V");  // meridional component

                double W_transf2 = face->get_parameter("Ninja" + std::to_string(d2));   // transfert function
                double U_lib2 = face->get_parameter("Ninja" + std::to_string(d2) + "_U");  // zonal component
                double V_lib2 = face->get_parameter("Ninja" + std::to_string(d2) + "_V");  // meridional component

                // Determine wind component from the wind field library using a weighted mean
                U = U_lib1*(theta2-theta)/(theta2-theta1)+U_lib2*(theta-theta1)/(theta2-theta1);
                V = V_lib1*(theta2-theta)/(theta2-theta1)+V_lib2*(theta-theta1)/(theta2-theta1);
                W_transf = W_transf1*(theta2-theta)/(theta2-theta1)+W_transf2*(theta-theta1)/(theta2-theta1);
            }

            // NEW wind direction from the wind field library
            theta = math::gis::zonal2dir(U, V);

            // NEW wind intensity from the wind field library
            W_transf = std::min(W_transf,Max_spdup);         // Limit speed up value to avoid unrelistic values at crest top
            W_transf = std::max(W_transf,Min_spdup);         
            double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v) * W_transf;
            W = std::max(W, 0.1);
            face->set_face_data("Ninja_speed", W);

            // Update U and V wind components
            U = -W  * sin(theta);
            V = -W  * cos(theta);

            //go back from H_forc to reference
            W = Atmosphere::log_scale_wind(W,
                                           H_forc,  //  Reference height for GEM forcing and WindNinja wind field library
                                           Atmosphere::Z_U_R,  // UR is at our reference height
                                           0); // no canopy, no snow, but uses a snow roughness

            face->set_face_data("U_R", W);
            face->set_face_data("vw_dir", theta * 180.0 / M_PI);

            face->set_face_data("Ninja_u", U); // these are still H_forc
            face->set_face_data("Ninja_v", V);

            Vector_2 v_corr = math::gis::bearing_to_cartesian(theta * 180.0 / M_PI);
            Vector_3 v3(-v_corr.x(), -v_corr.y(), 0); //negate as direction it's blowing instead of where it is from!!
            face->set_face_vector("wind_direction", v3);

            Vector_2 v_orig = math::gis::bearing_to_cartesian(theta_orig* 180.0 / M_PI);
            Vector_3 v3_orig(-v_orig.x(),-v_orig.y(), 0); //negate as direction it's blowing instead of where it is from!!
            face->set_face_vector("wind_direction_original",v3_orig);

            face->set_face_data("vw_dir_orig", theta_orig * 180.0 / M_PI);

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
                    u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(), neigh->face_data("U_R")));
            }

            auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());

            double new_u = face->get_module_data<data>(ID)->interp_smoothing(u, query);
            face->get_module_data<data>(ID)->temp_u = new_u;
        }
        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);
            face->set_face_data("U_R", std::max(0.1, face->get_module_data<data>(ID)->temp_u));
        }
   }

Ninja_wind::~Ninja_wind()
{

}
