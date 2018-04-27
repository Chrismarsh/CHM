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

#include "Kunkel_monthlyTd_rh.hpp"


Kunkel_monthlyTd_rh::Kunkel_monthlyTd_rh(config_file cfg)
        :module_base(parallel::data)

{
    provides("rh");
    provides("Td_lapse_rate");

    depends("t");
    depends_from_met("rh");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

Kunkel_monthlyTd_rh::~Kunkel_monthlyTd_rh()
{

}
void Kunkel_monthlyTd_rh::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void Kunkel_monthlyTd_rh::run(mesh_elem& face)
{
//    size_t ID = face->_debug_ID;
    // 1/km
    double lapse_rates[] = {
            0.41,
            0.42,
            0.40,
            0.39,
            0.38,
            0.36,
            0.33,
            0.33,
            0.36,
            0.37,
            0.4,
            0.4
    } ;

    double lapse = lapse_rates[ global_param->month() - 1 ] / 1000.; // -> 1/m

    //taken from mio
    const double  Bw = 17.502, Cw = 240.97; //parameters for water
    const double Bi = 22.452, Ci = 272.55; //parameters for ice

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("t")) || is_nan(s->get("rh")))
            continue;

        double t = s->get("t")+273.15;
        double rh = s->get("rh")/100.;

        double Tdz0 = mio::Atmosphere::RhtoDewPoint(rh,t,false) - 273.15; // K

        double C = t < 273.15 ? Ci : Cw;
        double B = t < 273.15 ? Bi : Bw;

        double z = 0.;
        double z0 = face->get_z();
//        double am = lapse;
        double Td_z = -lapse*C*(z-z0) / B + Tdz0;
//        double Td_z = (-am*(z-z0)*(C+Tdz0)/B+Tdz0)/(1+am*(z-z0)*(C+Tdz0)/(B*C));
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), Td_z ) );
    }

   

    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double Tdz0 = face->get_module_data<data>(ID)->interp(lowered_values, query);//C

    //raise value back up to the face's elevation from sea level
    double t = face->face_data("t") + 273.15;
    double C = t < 273.15 ? Ci : Cw;
    double B = t < 273.15 ? Bi : Bw;

    double z0 = 0.;
    double z = face->get_z();
//    double am = lapse;
    double Td_z = -lapse*C*(z-z0) / B + Tdz0;
//    double Td_z = (-am*(z-z0)*(C+Tdz0)/B+Tdz0)/(1+am*(z-z0)*(C+Tdz0)/(B*C));

    double rh = mio::Atmosphere::DewPointtoRh(Td_z+273.15,t,false);

    face->set_face_data("rh",rh*100.0);
    face->set_face_data("Td_lapse_rate",lapse);

}