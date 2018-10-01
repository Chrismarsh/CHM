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

#include "Dodson_NSA_ta.h"
REGISTER_MODULE_CPP(Dodson_NSA_ta);

Dodson_NSA_ta::Dodson_NSA_ta(config_file cfg)
:module_base(parallel::data)
{

    depends_from_met("t");
    provides("t");
    provides("t_lapse_rate");
}
Dodson_NSA_ta::~Dodson_NSA_ta()
{

}
void Dodson_NSA_ta::init(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);
        d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
    }
}
void Dodson_NSA_ta::run(mesh_elem &face)
{

    double Po = 100000.0; //sea level pressure (Pa)
    double Cp = 1005.0; //specific heat of dry air J/[KgK]
    double R = 8.3143; //gas constant J/(molK)
    double m = 0.02897; //molecular weight of dry air kg/mol
    double Tb = 300.0; //assumed sea level temperature 300K
    double g = 9.810616; //gravity m/s/s
    double lapse = 0.0065; //K/m

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;


    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("t")))
            continue;

        double elev = s->z(); //station elevation
        //pressure at station's elevation
        double Pz = Po * pow(Tb/(Tb+lapse*elev),(m*g)/(lapse*R));
        double ta = s->get("t") + 273.15; //to K

        //calculate virtual temp (eqn 1)
        double ratio = (Po/Pz);
        double exp = R/(m*Cp);
        double theta = ta * pow(ratio,exp);

        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), theta ) );
    }

    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());

    //interpolated virtual temp, now go back to station
    double theta = face->get_module_data<data>(ID)->interp(lowered_values, query);
    double elev = face->get_z();
    double Pz = Po * pow(Tb/(Tb+ (-lapse)*elev),(m*g)/((-lapse)*R));
    double ratio = (Po/Pz);
    double exp = R/(m*Cp);

    double Ta = ( theta/pow(ratio,exp) );
    Ta -= 273.15;

    face->set_face_data("t",Ta);
    face->set_face_data("t_lapse_rate",lapse);
}
