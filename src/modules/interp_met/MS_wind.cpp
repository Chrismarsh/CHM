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

#include "MS_wind.hpp"
REGISTER_MODULE_CPP(MS_wind);

MS_wind::MS_wind(config_file cfg)
        : module_base("MS_wind", parallel::domain, cfg)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("W_speedup");
    provides("vw_dir");
    provides("2m_zonal_u");
    provides("2m_zonal_v");

    provides("interp_zonal_u");
    provides("interp_zonal_v");

    provides("lookup_d");

    provides("vw_dir_orig");

//    provides_vector("wind_direction");

    speedup_height = cfg.get("speedup_height",2.0);
    use_ryan_dir = cfg.get("use_ryan_dir",false);
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void MS_wind::init(mesh& domain)
{
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

         auto d = face->make_module_data<data>(ID);
         d->interp.init(global_param->interp_algorithm,face->stations().size() );
         d->interp_smoothing.init(interp_alg::tpspline,3,{ {"reuse_LU","true"}});
    }
}


void MS_wind::run(mesh& domain)
{
    if(!use_ryan_dir)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);
		     std::vector<boost::tuple<double, double, double> > u;
		     std::vector<boost::tuple<double, double, double> > v;
		     for (auto &s : face->stations())
		     {
		       if (is_nan((*s)["U_R"_s]) || is_nan((*s)["vw_dir"_s]))
			 continue;

		       double theta = (*s)["vw_dir"_s] * M_PI / 180.;

		       double W = (*s)["U_R"_s];
		       W = std::max(W, 0.1);

		       W = Atmosphere::log_scale_wind(W,
						      Atmosphere::Z_U_R,  // UR is at our reference height
						      speedup_height,  // MS assumes a 2m wind speed
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

		     (*face)["interp_zonal_u"_s]= zonal_u;
		     (*face)["interp_zonal_v"_s]= zonal_v;

		     //Get back the interpolated wind direction
		     // -- not sure if there is a better way to do this, but at least a first order to getting the right direction
		     double theta = math::gis::zonal2dir(zonal_u, zonal_v);

		     double theta_orig = theta;

		     // Use this wind dir to figure out which lookUP we need
		     int d = int( std::round(theta * 180.0 / M_PI / 45.));

		     if (d == 8) d = 0; // floor(360/45) = 8, which we don't have, as 0 is already North, so use that.

		     (*face)["lookup_d"_s]= d;

		     // get the speedup for the interpolated direction
		     double U_speedup = face->parameter("MS" + std::to_string(d) + "_U");
		     double V_speedup = face->parameter("MS" + std::to_string(d) + "_V");
		     double W_speedup = face->parameter("MS" + std::to_string(d));

		     // Speed up interpolated zonal_u & zonal_v
		     double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v) * W_speedup;
		     W = std::max(W, 0.1);
		     (*face)["W_speedup"_s]= W;

		     //Now recover a U and V from this to get a new direction
		     double U = U_speedup * W;
		     double V = V_speedup * W;

		     // NEW wind direction after we do the U,V speedup
		     theta = math::gis::zonal2dir(U, V);

		     //go back from 2m to reference
		     W = Atmosphere::log_scale_wind(W,
						    speedup_height,  // MS assumes a 2m wind speed
						    Atmosphere::Z_U_R,  // UR is at our reference height
						    0); // no canopy, no snow, but uses a snow roughness

		     (*face)["U_R"_s]= W;
		     (*face)["vw_dir"_s]= theta * 180.0 / M_PI;

		     (*face)["2m_zonal_u"_s]= U; // these are still 2m
		     (*face)["2m_zonal_v"_s]= V;

		     Vector_2 v_corr = math::gis::bearing_to_cartesian(theta * 180.0 / M_PI);
		     Vector_3 v3(-v_corr.x(), -v_corr.y(), 0); //negate as direction it's blowing instead of where it is from!!
//		     face->set_face_vector("wind_direction", v3);


		     //        Vector_2 v_orig = math::gis::bearing_to_cartesian(theta_orig* 180.0 / M_PI);
		     //        Vector_3 v3_orig(-v_orig.x(),-v_orig.y(), 0); //negate as direction it's blowing instead of where it is from!!
		     //        face->set_face_vector("wind_direction_original",v3_orig);

		     (*face)["vw_dir_orig"_s]= theta_orig * 180.0 / M_PI;

        }

	// Need to access U_R from neighbors
	domain->ghost_neighbors_communicate_variable("U_R"_s);

        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
          auto face = domain->face(i);

          std::vector<boost::tuple<double, double, double>> u;
          for (size_t j = 0; j < 3; j++)
          {
            auto neigh = face->neighbor(j);

            if (neigh != nullptr)
            {
                try
                {
                    u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(), (*neigh)["U_R"_s]));
                }
                catch (...)
                {
                    LOG_DEBUG << "face global id " << face->cell_global_id;
                    LOG_DEBUG << "face is ghost? " << face->is_ghost;
                    LOG_DEBUG << "neigh that caused problem has global id " << neigh->cell_global_id;
                    LOG_DEBUG << "face neighbors are " << face->is_ghost;
                    for(int i = 0; i < 3; i++)
                    {
                        auto neigh = face->neighbor(i);
                        if (neigh != nullptr)
                        {
                            LOG_DEBUG << "\tneigh " << i << " global id " << neigh->cell_global_id;
                            LOG_DEBUG << "\tneigh " << i << " local id " << neigh->cell_local_id;
                            LOG_DEBUG << "\tneigh " << i << " is ghost? " << neigh->is_ghost;
                        }
                        else
                        {
                            LOG_DEBUG << "\tneigh " << i << " is nullptr";
                        }
                    }

                    CHM_THROW_EXCEPTION(module_error, "RIP");
                }
            }

          }

          double new_u = (*face)["U_R"_s];

          if (u.size() > 0)
          {
            auto query =
                boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
            new_u = face->get_module_data<data>(ID)->interp_smoothing(u, query);
          }

          face->get_module_data<data>(ID)->temp_u = new_u;
        }



        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);

		     (*face)["U_R"_s]= std::max(0.1, face->get_module_data<data>(ID)->temp_u);
        }

    }else
    {
        // omega_s needs to be scaled on [-0.5,0.5]
        double max_omega_s = -99999.0;

        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);

             std::vector<boost::tuple<double, double, double> > u;
             std::vector<boost::tuple<double, double, double> > v;
             for (auto &s : face->stations())
             {
               if (is_nan((*s)["U_R"_s]) || is_nan((*s)["vw_dir"_s]))
                 continue;

               double theta = (*s)["vw_dir"_s] * M_PI / 180.;

               auto f = domain->find_closest_face(s->x(),s->y());
               //figure out which lookup map we need
               int d = int(theta*180/M_PI/45.);
               if (d == 0) d = 8;
               double speedup = f->parameter("MS"+std::to_string(d));

               double W = (*s)["U_R"_s] / speedup;
               W = std::max(W, 0.1);
               W = Atmosphere::log_scale_wind(W,
                                              Atmosphere::Z_U_R,  // UR is at our reference height
                                              speedup_height,  // MS assumes a 2m wind speed
                                              0); // no canopy, no snow, but uses a snow roughness


               double zonal_u = -W * sin(theta);
               double zonal_v = -W * cos(theta);

               u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
               v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
             }
             //http://mst.nerc.ac.uk/wind_vect_convs.html

             auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
             double zonal_u = face->get_module_data<data>(ID)->interp(u, query);
             double zonal_v = face->get_module_data<data>(ID)->interp(v, query);

             double theta = 3.0 * M_PI * 0.5 - atan2(zonal_v, zonal_u);

             if (theta > 2.0 * M_PI)
               theta = theta - 2.0 * M_PI;

             //eqn 15
             double omega_s = face->slope() * cos(theta - face->aspect());

             if (fabs(omega_s) > max_omega_s)
               max_omega_s = fabs(omega_s);

             double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);

             face->get_module_data<data>(ID)->corrected_theta = theta;
             face->get_module_data<data>(ID)->W = W;

        }


        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);


		     double theta= face->get_module_data<data>(ID)->corrected_theta;
		     double W= face->get_module_data<data>(ID)->W;

		     //what liston calls 'wind slope'
		     double omega_s = face->slope() * cos(theta - face->aspect());

		     //scale between [-0.5,0.5]
		     omega_s = omega_s / (max_omega_s * 2.0);

		     //scale wind here via MS

		     double aspect = face->aspect();
		     double dirdiff = 0;

		     double d2r = M_PI/180.0;


		     //Ryan terain divergence
		     // follow LIston's micromet implimentation to handle 0-360
		     if (aspect > 270.0*d2r &&  theta < 90.0*d2r)
		       dirdiff = aspect - theta - 360.0*d2r;
		     else if (aspect < 90.0*d2r && theta > 270.0*d2r)
		       dirdiff = aspect - theta + 360.0*d2r;
		     else
		       dirdiff = aspect - theta;

		     if (std::abs(dirdiff) < 90.0*d2r)
		     {
		       theta = theta - 0.5 * std::min(omega_s, 45.0*d2r) * sin((2.0 * dirdiff));
		       if (theta > 2.0*M_PI)
			 theta = theta - 2*M_PI;
		       else if (theta < 0.0)
			 theta = theta + 2.0*M_PI;

		     }

		     //figure out which lookup map we need
		     int d = int(theta*180.0/M_PI/45.);
		     if (d == 0) d = 8;

		     double speedup = face->parameter("MS"+std::to_string(d));
		     W = W*speedup;

		     W = std::max(W,0.1);
		     W = std::min(W,30.0);

		     W = Atmosphere::log_scale_wind(W,
						    speedup_height,  // MS assumes a 2m wind speed
						    Atmosphere::Z_U_R,  // UR is at our reference height
						    0); // no canopy, no snow, but uses a snow roughness


		     (*face)["U_R"_s]= W;
		     (*face)["vw_dir"_s]= theta * 180.0 / M_PI;


		     Vector_2 v = math::gis::bearing_to_cartesian(theta* 180.0 / M_PI);
		     Vector_3 v3(-v.x(),-v.y(), 0); //negate as direction it's blowing instead of where it is from!!

//		     face->set_face_vector("wind_direction",v3);

        }

    // Need to access U_R from neighbors
    domain->ghost_neighbors_communicate_variable("U_R"_s);

    #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {

            auto face = domain->face(i);

		     std::vector<boost::tuple<double, double, double> > u;
		     for (size_t j = 0; j < 3; j++)
		     {
		       auto neigh = face->neighbor(j);
		       if (neigh != nullptr)
			 u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(),(*neigh)["U_R"_s]));
		     }


		     double new_u = (*face)["U_R"_s];
		     if (u.size() > 0)
		     {
		       auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
		       new_u = face->get_module_data<data>(ID)->interp_smoothing(u, query);
		     }

		     face->get_module_data<data>(ID)->temp_u = new_u;

        }

    #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);

		     (*face)["U_R"_s]= std::max(0.1,face->get_module_data<data>(ID)->temp_u) ;

        }

    }
}

//Old version, uses ryan. Deal with enabling this later, but for now use the version that gets dir off the u,v components of MS

//void MS_wind::run(mesh& domain)
//{
//
//    // omega_s needs to be scaled on [-0.5,0.5]
//    double max_omega_s = -99999.0;
//
//    #pragma omp parallel for
//    for (size_t i = 0; i < domain->size_faces(); i++)
//    {
//        auto face = domain->face(i);
//
//        std::vector<boost::tuple<double, double, double> > u;
//        std::vector<boost::tuple<double, double, double> > v;
//        for (auto &s : face->stations())
//        {
//            if (is_nan((*s)["U_R"_s]) || is_nan((*s)["vw_dir"_s]))
//                continue;
//
//            double theta = (*s)["vw_dir"_s] * M_PI / 180.;
//
//            auto f = domain->face(s->closest_face());
//            //figure out which lookup map we need
//            int d = int(theta*180/M_PI/45.);
//            if (d == 0) d = 8;
//            double speedup = f->parameter("MS"+std::to_string(d));
//
//            double W = (*s)["U_R"_s] / speedup;
//            W = std::max(W, 0.1);
//
//
//            double zonal_u = -W * sin(theta);
//            double zonal_v = -W * cos(theta);
//
//            u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
//            v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
//        }
////http://mst.nerc.ac.uk/wind_vect_convs.html
//
//        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
//        double zonal_u = face->get_module_data<data>(ID)->interp(u, query);
//        double zonal_v = face->get_module_data<data>(ID)->interp(v, query);
//
//        double theta = 3.0 * M_PI * 0.5 - atan2(zonal_v, zonal_u);
//
//        if (theta > 2.0 * M_PI)
//            theta = theta - 2.0 * M_PI;
//
//        //eqn 15
//        double omega_s = face->slope() * cos(theta - face->aspect());
//
//        if (fabs(omega_s) > max_omega_s)
//            max_omega_s = fabs(omega_s);
//
//        double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);
//
//        face->get_module_data<data>(ID)->corrected_theta = theta;
//        face->get_module_data<data>(ID)->W = W;
//    }
//    #pragma omp parallel for
//    for (size_t i = 0; i < domain->size_faces(); i++)
//    {
//        auto face = domain->face(i);
//
//        double theta= face->get_module_data<data>(ID)->corrected_theta;
//        double W= face->get_module_data<data>(ID)->W;
//
//        //what liston calls 'wind slope'
//        double omega_s = face->slope() * cos(theta - face->aspect());
//
//        //scale between [-0.5,0.5]
//        omega_s = omega_s / (max_omega_s * 2.0);
//
//        //scale wind here via MS
//
//        double aspect = face->aspect();
//        double dirdiff = 0;
//
//        double d2r = M_PI/180.0;
//
//
//        //Ryan terain divergence
//        // follow LIston's micromet implimentation to handle 0-360
//        if (aspect > 270.0*d2r &&  theta < 90.0*d2r)
//            dirdiff = aspect - theta - 360.0*d2r;
//        else if (aspect < 90.0*d2r && theta > 270.0*d2r)
//            dirdiff = aspect - theta + 360.0*d2r;
//        else
//            dirdiff = aspect - theta;
//
//        if (std::abs(dirdiff) < 90.0*d2r)
//        {
//            theta = theta - 0.5 * std::min(omega_s, 45.0*d2r) * sin((2.0 * dirdiff));
//            if (theta > 2.0*M_PI)
//                theta = theta - 2*M_PI;
//            else if (theta < 0.0)
//                theta = theta + 2.0*M_PI;
//
//        }
//
//        //figure out which lookup map we need
//        int d = int(theta*180.0/M_PI/45.);
//        if (d == 0) d = 8;
//
//        double speedup = face->parameter("MS"+std::to_string(d));
//        W = W*speedup;
//
//        W = std::max(W,0.1);
//        W = std::min(W,30.0);
//        (*face)["U_R"_s]= W;
//        (*face)["vw_dir"_s]= theta * 180.0 / M_PI;
//
//
//        Vector_2 v = math::gis::bearing_to_cartesian(theta* 180.0 / M_PI);
//        Vector_3 v3(-v.x(),-v.y(), 0); //negate as direction it's blowing instead of where it is from!!
//
//        face->set_face_vector("wind_direction",v3);
//
//    }
//
//
//#pragma omp parallel for
//    for (size_t i = 0; i < domain->size_faces(); i++)
//    {
//
//        auto face = domain->face(i);
//        std::vector<boost::tuple<double, double, double> > u;
//        for (size_t j = 0; j < 3; j++)
//        {
//            auto neigh = face->neighbor(j);
//            if (neigh != nullptr)
//                u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(),(*neigh)["U_R"_s]));
//        }
//
//        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
//
//        double new_u = face->get_module_data<data>(ID)->interp_smoothing(u, query);
//        face->get_module_data<data>(ID)->temp_u = new_u;
//    }
//#pragma omp parallel for
//    for (size_t i = 0; i < domain->size_faces(); i++)
//    {
//        auto face = domain->face(i);
//        (*face)["U_R"_s]= std::max(0.1,face->get_module_data<data>(ID->temp_u) );
//    }
//}

MS_wind::~MS_wind()
{

}
