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

#include "Thornton_p.hpp"

Thornton_p::Thornton_p(config_file cfg)
        :module_base(parallel::data)
{

    depends_from_met("p");

    provides("p");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}
void Thornton_p::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void Thornton_p::run(mesh_elem& face)
{
    //km^-1
    double monthly_factors[] = {0.35, 0.35, 0.35, 0.30, 0.25, 0.20, 0.20, 0.20, 0.20, 0.25, 0.30, 0.35};
    for(auto& mf: monthly_factors)
    {
        mf /= 1000.0; //to m^-1
    }
    std::vector< boost::tuple<double, double, double> > ppt;
    std::vector< boost::tuple<double, double, double> > staion_z;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("p")))
            continue;
        double u = s->get("p");
        ppt.push_back( boost::make_tuple(s->x(), s->y(), u ) );
        staion_z.push_back( boost::make_tuple(s->x(), s->y(), s->z() ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double p0 = face->get_module_data<data>(ID)->interp(ppt, query);
    double z0 = face->get_module_data<data>(ID)->interp(staion_z,query);
    double z = face->get_z();

    int month = global_param->month()-1;
    double lapse = monthly_factors[month];

    double P = p0*( (1+lapse*(z-z0))/(1-lapse*(z-z0)));
    P = std::max(0.0,P);

    face->set_face_data("p", P);


}

Thornton_p::~Thornton_p()
{

}