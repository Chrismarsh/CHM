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

#include "Dist_tlapse.hpp"
REGISTER_MODULE_CPP(Dist_tlapse);

Dist_tlapse::Dist_tlapse(config_file cfg)
        : module_base("Dist_tlapse", parallel::data, cfg)

{
    provides("t");
    provides("t_lapse_rate");

    depends_from_met("t");
    depends_from_met("t_lapse_rate");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Dist_tlapse::~Dist_tlapse()
{

}
void Dist_tlapse::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void Dist_tlapse::run(mesh_elem& face)
{
    // Distributed forcing lapse_rate to each face (changes per time step)
    //lapse_rate=0.0047;

    // Find nearest station
    auto s_near = global_param->nearest_station(face->center().x(),face->center().y());

    // Get the lapse rate from that station
    double lapse_rate = s_near.at(0)->get("t_lapse_rate");

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("t")))
            continue;
        double v = s->get("t") - lapse_rate * (0.0 - s->z());
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double value = face->get_module_data<data>(ID)->interp(lowered_values, query);

    //raise value back up to the face's elevation from sea level
    value =  value + lapse_rate * (0.0 - face->get_z());

    face->set_face_data("t",value);

    face->set_face_data("t_lapse_rate",lapse_rate);

}
