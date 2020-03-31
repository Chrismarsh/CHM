#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
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

#include "WindNinja.hpp"
REGISTER_MODULE_CPP(WindNinja);

WindNinja::WindNinja(config_file cfg)
        : module_base("WindNinja", parallel::domain, cfg)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("Ninja_speed");
    provides("Ninja_speed_nodown");

    provides("vw_dir");
    provides("Ninja_u");
    provides("Ninja_v");

    provides("omega_s");
    provides("W_transf");

    provides("interp_zonal_u");
    provides("interp_zonal_v");

    provides("lookup_d");

    provides("vw_dir_orig");

    provides_vector("wind_direction_original");
    provides_vector("wind_direction");

    ninja_average = cfg.get("ninja_average",true);

    compute_Sx = cfg.get("compute_Sx",true);
    // We are going to use the winstral parameterization of Sx to modify out windfield. So we cannot do it later
    // it needs to be coupled with  the windspeed parameterization. Manually instantiate our own copy of this module
    // and we can then access Sx as if it were any other object. The config for this will need to come from WN's config section
    if(compute_Sx)
    {
        optional("snowdepthavg");
        provides("Sx");
        conflicts("Winstral_parameters"); // we cannot have Winstral_parameters alongside this if we generate the Sx.

        config_file tmp;
        tmp.put("angular_window",30.);
        tmp.put("size_of_step",10.);
        Sx = boost::dynamic_pointer_cast<Winstral_parameters>(module_factory::create("Winstral_parameters",tmp));
    }

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void WindNinja::init(mesh& domain)
{
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,face->stations().size() );
        d->interp_smoothing.init(interp_alg::tpspline,3,{ {"reuse_LU","true"}});
    }

    N_windfield = 0;
    for(auto& itr: domain->parameters() )
    {
        if( itr.find("Ninja") != std::string::npos &&
            itr.find("_U") != std::string::npos) // ignore all _U, _V, and _LAvg params
        {
            ++N_windfield;
        }
    }

    LOG_DEBUG << "Found " << N_windfield << " windfields";
    try
    {
        // see if we have a manually specified number, error out as this is deprecated
        // favour error as we dont want the user to think something is working when it's not
        cfg.get<int>("N_windfield");
        CHM_THROW_EXCEPTION(module_error, "N_windfield is deprecated and no longer used");
    }
    catch(pt::ptree_bad_path& e)
    {
        // not user specified, is ok
    }

    if (N_windfield == 0)
    {
        CHM_THROW_EXCEPTION(module_error,"Could not find any required wind ninja maps");
    }

    // ensure a user-specified L_avg value makes sense
    L_avg = -1;
    try
    {
        L_avg = cfg.get<int>("L_avg");

        int n_lavg = 0; // number of params that have a specific l_avg
        auto s_lavg = std::to_string(L_avg);

        for(auto& itr: domain->parameters() )
        {
            if( itr.find("Ninja") != std::string::npos &&
                 itr.find("_"+s_lavg) != std::string::npos
            )
            {
                ++n_lavg;
            }
        }

        if(n_lavg != N_windfield)
        {
            CHM_THROW_EXCEPTION(module_error,"The user specified value of L_avg=" +
                s_lavg + " does not match the number of wind fields.");
        }
    }
    catch(...)
    {
        // not user specified, is ok
    }

    // We need to handle the following case:
    // if L_avg is not defined in the cfg file BUT the mesh contains only WN parameters with an explicit values for L_avg.
    if(L_avg == -1)
    {
        std::set<int> param_found_L_avg={};
        for(auto& itr: domain->parameters() )
        {
            if( itr.find("Ninja") != std::string::npos &&
                itr.find("_") != std::string::npos  &&
                itr.find("_U") == std::string::npos && // AND these are the U and V fields as the Lavg isn't on U or V
                itr.find("_V") == std::string::npos
                )
            {
                // see if we can auto detect an Lavg value
                auto pos = itr.find("_") + 1;
                auto param_L_avg = -1;
                try
                {
                    param_L_avg = std::stoi(itr.substr(pos));
                    param_found_L_avg.insert(param_L_avg);
                }
                catch(std::invalid_argument& e)
                {

                    CHM_THROW_EXCEPTION(module_error,"Found a WindNinja L_avg parameter that has a malformed name. Param = " + itr + ". L_avg = " + itr.substr(pos) );
                }


            }
        }

        if(param_found_L_avg.size() == 1) // we found only 1 possible Lavg set of values, so we can use it.
        {
            L_avg = *(param_found_L_avg.begin());
            LOG_DEBUG << "Using L_avg from parameters. L_avg = " << L_avg;
        }
        else
        {
            std::string found = "[   ";
            for(auto& f : param_found_L_avg)
            {
                found += std::to_string(f);
                found += "   ";
            }
            found += " ]";
            // Ok we have a problem. There are multiple NinjaX_LAVG params in the file
            CHM_THROW_EXCEPTION(module_error,"WindNinja: Multiple parameters with possible L_avg values found, but L_avg was not specified. "
                                             "Please set L_avg in the config file. Found the following L_avg values: " + found );

        }
    }

    // Ensure that L_avg, either given or found, is valid
    int valid_Lavg_params=0;

    for(int theta = 0; theta<=2*M_PI;theta++)
    {
        auto delta_angle = 360. / N_windfield;
        auto d = int(theta * 180.0 / M_PI / delta_angle);
        if (d == 0) d = N_windfield;
        auto face = domain->face(0);

        std::string param="";

        if(L_avg == -1)
        {
            param = "Ninja" + std::to_string(d);
        }else
        {
            param = "Ninja" + std::to_string(d) +'_' + std::to_string(L_avg);   // transfert function
        }

        if( !face->has_parameter(param))
        {
            CHM_THROW_EXCEPTION(module_error,"WindNinja: Missing parameter: " + param);
        }

        param = "Ninja" + std::to_string(d) + "_U";
        if( !face->has_parameter(param))
        {
            CHM_THROW_EXCEPTION(module_error,"WindNinja: Missing parameter: " + param);
        }

        param = "Ninja" + std::to_string(d) + "_V";
        if( !face->has_parameter(param))
        {
            CHM_THROW_EXCEPTION(module_error,"WindNinja: Missing parameter: " + param);
        }

    }


    H_forc = cfg.get("H_forc",40.0);
    Max_spdup = cfg.get("Max_spdup",3.);
    Min_spdup = cfg.get("Min_spdup",0.1);
    ninja_recirc = cfg.get("ninja_recirc",false);
    Sx_crit = cfg.get("Sx_crit", 30.);
}


void WindNinja::run(mesh& domain)
{
        double transf_max = -9999.0;
        double max_omega_s = -9999.0;

        double delta_angle = 360. / N_windfield;

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
                                               H_forc,  //  Reference height for GEM forcing and WindNinja wind field library
                                               0); // no canopy, no snow, but uses a snow roughness

                double zonal_u = -W * sin(theta);
                double zonal_v = -W * cos(theta);

                u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
                v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
            }

            /**
             * Example use of Sx, modify as needed in the future
             * When we are reusing other modules' members they ofc need to be thread safe
             */
//            if(compute_Sx)
//                (*face)["Sx"_s]= Sx->Sx(domain,face);

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

            // Save original direction
            Vector_2 v_orig = math::gis::bearing_to_cartesian(theta_orig* 180.0 / M_PI);
            Vector_3 v3_orig(-v_orig.x(),-v_orig.y(), 0); //negate as direction it's blowing instead of where it is from!!
            face->set_face_vector("wind_direction_original",v3_orig);
            (*face)["vw_dir_orig"_s]= theta_orig * 180.0 / M_PI;

            double U = 0.;
            double V = 0.;
            double W_transf = 0.;


            if(!ninja_average)  // No Linear interpolation between the closest 2 wind fields from the library
            {
                // Use this wind dir to figure out which wind field from the library we need
                // Wind field are available each delta_angle deg.
                int d = int(theta * 180.0 / M_PI / delta_angle);
                if (d == 0) d = N_windfield;
                (*face)["lookup_d"_s]= d;

                // get the transfert function and associated wind component for the interpolated wind direction
		if(L_avg == -1)
                {		                 
			W_transf = face->parameter("Ninja" + std::to_string(d));   // transfert function
		}else
		{
			W_transf = face->parameter("Ninja" + std::to_string(d) +'_' + std::to_string(L_avg));   // transfert function
                }
                U = face->parameter("Ninja" + std::to_string(d) + "_U");  // zonal component
                V = face->parameter("Ninja" + std::to_string(d) + "_V");  // meridional component

           }else // Linear interpolation between the closest 2 wind fields from the library
           {

                // Use this wind dir to figure out which wind fields from the library we need
                // Wind fields are available each 15 deg.
                int d1 = int(theta * 180.0 / M_PI / delta_angle);
                double theta1 = d1 * delta_angle * M_PI / 180.0;
                if (d1 == 0) d1 = N_windfield;

                int d2 = int((theta * 180.0 / M_PI + delta_angle) / delta_angle);
                double theta2 = d2 * delta_angle *M_PI / 180.0;
                if (d2 == 0) d2 = N_windfield;

                double d = d1*(theta2-theta)/(theta2-theta1)+d2*(theta-theta1)/(theta2-theta1);
                (*face)["lookup_d"_s]= d;

		double W_transf1 = 0.;
                // get the transfert function and associated wind component for the interpolated wind direction
		if(L_avg == -1)
		{	
                	W_transf1 = face->parameter("Ninja" + std::to_string(d1));   // transfert function
                }else
		{
			W_transf1 = face->parameter("Ninja" + std::to_string(d1)+'_' + std::to_string(L_avg));   // transfert function
		}

                double U_lib1 = face->parameter("Ninja" + std::to_string(d1) + "_U");  // zonal component
                double V_lib1 = face->parameter("Ninja" + std::to_string(d1) + "_V");  // meridional component

		double W_transf2 = 0.;
		if(L_avg == -1)
		{	
                	W_transf2 = face->parameter("Ninja" + std::to_string(d2));   // transfert function
                }else
		{
			W_transf2 = face->parameter("Ninja" + std::to_string(d2)+'_' + std::to_string(L_avg));   // transfert function
                }
		double U_lib2 = face->parameter("Ninja" + std::to_string(d2) + "_U");  // zonal component
                double V_lib2 = face->parameter("Ninja" + std::to_string(d2) + "_V");  // meridional component

                // Determine wind component from the wind field library using a weighted mean
                U = U_lib1*(theta2-theta)/(theta2-theta1)+U_lib2*(theta-theta1)/(theta2-theta1);
                V = V_lib1*(theta2-theta)/(theta2-theta1)+V_lib2*(theta-theta1)/(theta2-theta1);
                W_transf = W_transf1*(theta2-theta)/(theta2-theta1)+W_transf2*(theta-theta1)/(theta2-theta1);
            }

            if(fabs(W_transf)> transf_max )
                     transf_max  = fabs(W_transf);

            // NEW wind direction from the wind field library
            theta = math::gis::zonal2dir(U, V);

            //Compute what liston calls 'wind slope' using updated wind direction
     //       double omega_s = face->slope() * cos(theta - face->aspect());

     //       if (fabs(omega_s) > max_omega_s)
     //           max_omega_s = fabs(omega_s);

            double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);

            face->get_module_data<data>(ID)->corrected_theta = theta;
            face->get_module_data<data>(ID)->W = W;
            face->get_module_data<data>(ID)->W_transf = W_transf;

           }

        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
       {

           auto face = domain->face(i);

           double theta= face->get_module_data<data>(ID)->corrected_theta;
           double W= face->get_module_data<data>(ID)->W;
           double W_transf= face->get_module_data<data>(ID)->W_transf;

            (*face)["Ninja_speed_nodown"_s]= W;   // Wind speed without downscaling


           (*face)["vw_dir"_s]= theta * 180.0 / M_PI;
           // Limit speed up value to Max_spdup
           // Can be used to avoid unrelistic values at crest top
           if(W_transf>1. and transf_max>Max_spdup)
               W_transf = 1.+(Max_spdup-1.)*(W_transf-1.)/(transf_max-1.);

           if(compute_Sx)
           {
               if (ninja_recirc)
               {  // Need further test

              //Compute what liston calls 'wind slope' using updated wind direction
            //  double omega_s = face->slope() * cos(theta - face->aspect());
            //  (*face)["omega_s"_s]= omega_s;

         //     if( omega_s<-0.35 and W_transf>1.05)  //Reduce wind speed on the lee side of mountain crest
          //         W_transf = std::max(0.5,0.5+(omega_s+0.5)/(0.15)*0.5);   // Reduction 0f 50% for omega_s larger than 30 deg

                  double sx_loc = Sx->Sx(domain,face);
                  (*face)["Sx"_s] =sx_loc;
                  if( sx_loc>Sx_crit )  //Reduce wind speed on the lee side of mountain crest identified by Sx>Sx_crit
                       W_transf = 0.25;
                }
            }

            (*face)["W_transf"_s]= W_transf;

            // NEW wind intensity from the wind field library
            W = W * W_transf;
            W = std::max(W, 0.1);
            (*face)["Ninja_speed"_s]= W;    // Wind speed with downscaling

            // Update U and V wind components
            double U = -W  * sin(theta);
            double V = -W  * cos(theta);

            //go back from H_forc to reference
            W = Atmosphere::log_scale_wind(W,
                                           H_forc,  //  Reference height for GEM forcing and WindNinja wind field library
                                           Atmosphere::Z_U_R,  // UR is at our reference height
                                           0); // no canopy, no snow, but uses a snow roughness

            (*face)["U_R"_s]= W;


            (*face)["Ninja_u"_s]= U; // these are still H_forc
            (*face)["Ninja_v"_s]= V;

            Vector_2 v_corr = math::gis::bearing_to_cartesian(theta * 180.0 / M_PI);
            Vector_3 v3(-v_corr.x(), -v_corr.y(), 0); //negate as direction it's blowing instead of where it is from!!
            face->set_face_vector("wind_direction", v3);


        }


        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {

            auto face = domain->face(i);
            std::vector<boost::tuple<double, double, double> > u;
            for (size_t j = 0; j < 3; j++)
            {
                auto neigh = face->neighbor(j);
                if (neigh != nullptr && !neigh->_is_ghost)
                    u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(),(*neigh)["U_R"_s]));
            }

            auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
            if(u.size() > 0)
            {
                double new_u = face->get_module_data<data>(ID)->interp_smoothing(u, query);
                face->get_module_data<data>(ID)->temp_u = new_u;
            }
            else
            {
                face->get_module_data<data>(ID)->temp_u = (*face)["U_R"_s];
            }

        }
        #pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);
            (*face)["U_R"_s]= std::max(0.1, face->get_module_data<data>(ID)->temp_u);
        }
   }

WindNinja::~WindNinja()
{

}

#pragma clang diagnostic pop