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

#include "uniform_wind.hpp"
REGISTER_MODULE_CPP(uniform_wind);


uniform_wind::uniform_wind(config_file cfg)
        : module_base("uniform_wind", parallel::domain, cfg)

{
    depends_from_met("U_R");
    depends_from_met("vw_dir");

    provides("U_R");
    provides("vw_dir");

    provides_vector("wind_direction");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

//Calculates the curvature required
void uniform_wind::init(mesh& domain)
{

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<lwinddata>(ID);
        d.interp.init(global_param->interp_algorithm,face->stations().size() );
        face->coloured = false;
    }

}


void uniform_wind::run(mesh& domain)
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
           double zonal_u = -W * sin(theta);
           double zonal_v = -W * cos(theta);
           u.push_back(boost::make_tuple(s->x(), s->y(), zonal_u));
           v.push_back(boost::make_tuple(s->x(), s->y(), zonal_v));
        }

        // Interp over stations
        auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
        double zonal_u = face->get_module_data<lwinddata>(ID).interp(u, query);
        double zonal_v = face->get_module_data<lwinddata>(ID).interp(v, query);

        // Convert back to direction and magnitude
        double theta = 3.0 * M_PI * 0.5 - atan2(zonal_v, zonal_u);

        if (theta > 2 * M_PI)
         theta = theta - 2 * M_PI;

        double W = sqrt(zonal_u * zonal_u + zonal_v * zonal_v);
        double corrected_theta = 3.0 * M_PI * 0.5 - atan2(zonal_v, zonal_u);

        if (corrected_theta > 2 * M_PI)
         corrected_theta = corrected_theta - 2 * M_PI;

        face->get_module_data<lwinddata>(ID).corrected_theta = corrected_theta;
        face->get_module_data<lwinddata>(ID).W = W;

    }


    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        double corrected_theta= face->get_module_data<lwinddata>(ID).corrected_theta;
        double W= face->get_module_data<lwinddata>(ID).W;

        W = std::max(W,0.1);
        (*face)["U_R"_s]= W;
        (*face)["vw_dir"_s]= corrected_theta * 180.0 / M_PI;

        Vector_2 v_corr = math::gis::bearing_to_cartesian(corrected_theta * 180.0 / M_PI);
        Vector_3 v3(-v_corr.x(), -v_corr.y(), 0); //negate as direction it's blowing instead of where it is from!!
        face->set_face_vector("wind_direction", v3);
    }

}

uniform_wind::~uniform_wind()
{

}
