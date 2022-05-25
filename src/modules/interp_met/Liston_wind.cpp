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

#include "Liston_wind.hpp"
REGISTER_MODULE_CPP(Liston_wind);

Liston_wind::Liston_wind(config_file cfg)
        : module_base("Liston_wind", parallel::domain, cfg)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("U_R_orig");

    provides("zonal_u");
    provides("zonal_v");

    provides("vw_dir");
    provides("vw_dir_orig");
    provides("vw_dir_divergence");

    provides_vector("wind_direction");
    provides_vector("wind_direction_original");

    provides_parameter("Liston_curvature");
    distance = cfg.get<double>("distance",300);
    Ww_coeff = cfg.get<double>("Ww_coeff",1.0);

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void Liston_wind::init(mesh& domain)
{

    ys = cfg.get("ys",0.5);
    yc = cfg.get("yc",0.5);

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);
        auto& d = face->make_module_data<lwinddata>(ID);
        d.interp.init(global_param->interp_algorithm,face->stations().size() );
        d.interp_smoothing.init(interp_alg::tpspline,3,{ {"reuse_LU","true"}});

        face->coloured = false;

    }


    double curmax = -9999.0;

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);

        Point_3 me = face->center();
        mesh_elem north;
        mesh_elem south;
        mesh_elem east;
        mesh_elem west;

        mesh_elem northeast;
        mesh_elem northwest;
        mesh_elem southeast;
        mesh_elem southwest;

        north = domain->find_closest_face( math::gis::point_from_bearing(me,0,distance) ) ; // me.x(), me.y() + distance
        south = domain->find_closest_face( math::gis::point_from_bearing(me,180,distance) ); //me.x(), me.y() - distance
        west = domain->find_closest_face(  math::gis::point_from_bearing(me,270,distance) ); //me.x() - distance, me.y()
        east = domain->find_closest_face(  math::gis::point_from_bearing(me,90,distance)  ); //me.x() + distance, me.y()

        double z = face->get_z();
        double zw = west->get_z();
        double ze = east->get_z();
        double zs = south->get_z();
        double zn = north->get_z();

        double znw = 0.;
        double zne = 0.;
        double zse = 0.;
        double zsw = 0.;

        northeast = domain->find_closest_face(math::gis::point_from_bearing(me,45,distance)); //me.x() + distance, me.y() + distance
        zne = northeast->get_z();

        northwest = domain->find_closest_face(math::gis::point_from_bearing(me,315,distance)); //me.x() - distance, me.y() + distance
        znw = northwest->get_z();

        southeast = domain->find_closest_face(math::gis::point_from_bearing(me,135,distance)); //me.x() + distance, me.y() - distance
        zse = southeast->get_z();

        southwest = domain->find_closest_face(math::gis::point_from_bearing(me,225,distance)); //me.x() - distance, me.y() - distance
        zsw = southwest->get_z();

        double curve = .25 * ((z - .5 * (zw + ze)) / (2.0 * distance) + (z - .5 * (zs + zn)) / (2.0 * distance) +
                             (z - .5 * (zsw + zne)) / (2.0 * sqrt(2.0) * distance) +
                             (z - .5 * (znw + zse)) / (2.0 * sqrt(2.0) * distance));

        auto& c = face->get_module_data<lwinddata>(ID);
        c.curvature = curve;

        if (fabs(curve) > curmax)
        {
           curmax = fabs(curve);
        }

    }


    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);
        auto& c = face->get_module_data<lwinddata>(ID);

        double value = c.curvature / (curmax * 2.0);//rescale to [-0.5,+0.5];

        //with very coarse meshes, with very few total triangles,
        // there are edge cases where curmax=0 and makes curvature NAN. Just set it to 0, no curvature, and don't do silly speedup/down
        c.curvature = std::isnan(value) ? 0 : value;

        face->parameter("Liston_curvature"_s) =  c.curvature;

    }



//    if ( cfg.get("serialize",false) )
//    {
//        LOG_DEBUG << "Serializing liston curvature";
//        domain->serialize_parameter(cfg.get("serialize_output", "liston_curvature.mesh"),
//                                    "Liston curvature");
//    }

}


void Liston_wind::run(mesh& domain)
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

           double W = (*s)["U_R"_s];
           W = std::max(W, 0.1);

           double theta = (*s)["vw_dir"_s] * M_PI / 180.;
           double phi = math::gis::bearing_to_polar((*s)["vw_dir"_s] );

           double zonal_u = -W * sin(theta);//negate as it needs to be the direction the wind is *going*
           double zonal_v = -W * cos(theta);

           u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
           v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
        }
        //http://mst.nerc.ac.uk/wind_vect_convs.html

        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
        double zonal_u = face->get_module_data<lwinddata>(ID).interp(u, query);
        double zonal_v = face->get_module_data<lwinddata>(ID).interp(v, query);

        double theta = 3.0 * M_PI * 0.5 - atan2(zonal_v, zonal_u);
        //        double theta = atan2(-zonal_v, -zonal_u);
        if (theta > 2.0 * M_PI)
         theta = theta - 2.0 * M_PI;

        //eqn 15
        double omega_s = face->slope() * cos(theta - face->aspect());

        if (fabs(omega_s) > max_omega_s)
         max_omega_s = fabs(omega_s);

        double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);

        face->get_module_data<lwinddata>(ID).corrected_theta = theta;
        face->get_module_data<lwinddata>(ID).W = W;

        (*face)["U_R_orig"_s] = W;

        // Write updated U and V wind components
        double U = -W  * sin(theta);
        double V = -W  * cos(theta);
        (*face)["zonal_u"_s] = U;
        (*face)["zonal_v"_s] = V;

        // Save original direction
        Vector_2 v_orig = math::gis::bearing_to_cartesian(theta* 180.0 / M_PI);
        Vector_3 v3_orig(-v_orig.x(),-v_orig.y(), 0); //negate as direction it's blowing instead of where it is from!!
        face->set_face_vector("wind_direction_original",v3_orig);
        (*face)["vw_dir_orig"_s]= theta * 180.0 / M_PI;

    }


    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);

        double theta = face->get_module_data<lwinddata>(ID).corrected_theta;
        double W = face->get_module_data<lwinddata>(ID).W;

        //what liston calls 'wind slope'
        double omega_s = face->slope() * cos(theta - face->aspect());

        //scale between [-0.5,0.5]
        omega_s = omega_s / (max_omega_s * 2.0);

        double omega_c = face->parameter("Liston_curvature"_s);

        double Ww = Ww_coeff + ys * omega_s + yc * omega_c;

        if(std::isnan(Ww))
         Ww=Ww_coeff;

        W = W * Ww;


        //Ryan (1977) Terrain divergence
        // follow Liston's MicroMet implementation to handle 0-360
        // https://github.com/wk1984/Snow_Model_Fortran_Glen/blob/master/micromet_code.f#L1365

        double aspect = face->aspect();
        double dirdiff = 0;

        double d2r = M_PI/180.0;

        if (aspect > 270.0*d2r &&  theta < 90.0*d2r)
            dirdiff = aspect - theta - 360.0*d2r;
        else if (aspect < 90.0*d2r && theta > 270.0*d2r)
            dirdiff = aspect - theta + 360.0*d2r;
        else
            dirdiff = aspect - theta;

        if (std::abs(dirdiff) <= 90.0*d2r)
        {
           theta = theta - 0.5 * omega_s * sin((2.0 * dirdiff));
           if (theta > 2.0*M_PI)
             theta = theta - 2*M_PI;
           else if (theta < 0.0)
             theta = theta + 2.0*M_PI;
        }

        W = std::max(W,0.1);
        (*face)["U_R"_s]= W;
        (*face)["vw_dir"_s]= theta * 180.0 / M_PI;

        Vector_2 v = math::gis::bearing_to_cartesian(theta* 180.0 / M_PI);
        Vector_3 v3(-v.x(),-v.y(), 0); //negate as direction it's blowing instead of where it is from!!

        (*face)["vw_dir_divergence"_s]=fabs( -0.5 * omega_s * sin((2.0 * dirdiff)))*180.0/M_PI ;
        face->set_face_vector("wind_direction",v3);
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
             u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(), (*neigh)["U_R"_s]));
        }

        double new_u = (*face)["U_R"_s];
        if(u.size() > 0)
        {
           auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
           new_u = face->get_module_data<lwinddata>(ID).interp_smoothing(u, query);
        }

        face->get_module_data<lwinddata>(ID).temp_u = new_u;

    }


#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

        auto face = domain->face(i);
        (*face)["U_R"_s]=face->get_module_data<lwinddata>(ID).temp_u;

    }

}

Liston_wind::~Liston_wind()
{

}
