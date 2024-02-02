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

#include "t_no_lapse.hpp"
REGISTER_MODULE_CPP(t_no_lapse);

t_no_lapse::t_no_lapse(config_file cfg)
        : module_base("t_no_lapse", parallel::data, cfg)

{
    provides("t");
    provides("t_lapse_rate");

    depends_from_met("t");

    SPDLOG_DEBUG("Successfully instantiated module {}",this->ID);
}

t_no_lapse::~t_no_lapse()
{

}
void t_no_lapse::init(mesh& domain)
{

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<t_no_lapse::data>(ID);
        d.interp.init(global_param->interp_algorithm,face->stations().size() );
    }

}
void t_no_lapse::run(mesh_elem& face)
{

    // Const lapse rate
    double lapse_rate = 0.0;

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : face->stations())
    {
        if( is_nan((*s)["t"_s]))
            continue;
        double v = (*s)["t"_s] - lapse_rate * (0.0 - s->z());
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double value = face->get_module_data<data>(ID).interp(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    value =  value + lapse_rate * (0.0 - face->get_z());

    (*face)["t"_s]=value;
    (*face)["t_lapse_rate"_s]=lapse_rate;

}
