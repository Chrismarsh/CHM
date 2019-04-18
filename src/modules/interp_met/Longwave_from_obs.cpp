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

#include "Longwave_from_obs.hpp"
REGISTER_MODULE_CPP(Longwave_from_obs);

Longwave_from_obs::Longwave_from_obs(config_file cfg)
        : module_base("Longwave_from_obs", parallel::data, cfg)

{
    provides("ilwr");
    //provides("lw_lapse_rate");

    depends_from_met("Qli");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Longwave_from_obs::~Longwave_from_obs()
{

}
void Longwave_from_obs::init(mesh& domain)
{

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

	       auto face = domain->face(i);
	       auto d = face->make_module_data<data>(ID);
	       d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());

    }

}
void Longwave_from_obs::run(mesh_elem& face)
{

    // Use constant annual LW lapse rate based on emperical study in Alps
    double lapse_rate = 2.8/100; // 2.8 W/m^2 / 100 meters (Marty et al. 2002)

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("Qli")))
            continue;
        double v = s->get("Qli") - lapse_rate * (0.0 - s->z());
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double value = face->get_module_data<data>(ID)->interp(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    value =  value + lapse_rate * (0.0 - face->get_z());

    (*face)["ilwr"_s]=value;

    //(*face)["lw_lapse_rate"_s]=lapse_rate;

}
