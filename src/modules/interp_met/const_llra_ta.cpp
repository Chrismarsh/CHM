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

#include "const_llra_ta.hpp"
REGISTER_MODULE_CPP(const_llra_ta);

const_llra_ta::const_llra_ta(config_file cfg)
        : module_base("const_llra_ta", parallel::data, cfg)

{
    provides("t");
    provides("const_llra_ta");

    depends_from_met("t");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

const_llra_ta::~const_llra_ta()
{

}

void const_llra_ta::init(mesh& domain)
{

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

	       auto face = domain->face(i);
	       auto d = face->make_module_data<const_llra_ta::data>(ID);
	       d->interp.init(global_param->interp_algorithm,face->stations().size() );

    }
    LOG_DEBUG << "Successfully init module " << this->ID;

}
void const_llra_ta::run(mesh_elem& face)
{

    double lapse_rate = 0.0065;


    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : face->stations())
    {
        if( is_nan((*s)["t"]))
            continue;
        double v = (*s)["t"] - lapse_rate * (0.0 - s->z());
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double value = face->get_module_data<data>(ID)->interp(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    value =  value + lapse_rate * (0.0 - face->get_z());

    (*face)["t"_s]=value;
    (*face)["const_llra_ta"_s]=value;

}
