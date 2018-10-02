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

#include "Walcek_cloud.hpp"
REGISTER_MODULE_CPP(Walcek_cloud);

Walcek_cloud::Walcek_cloud(config_file cfg)
        : module_base("Walcek_cloud", parallel::data, cfg)
{
    provides("cloud_frac");
    depends("t");
    depends("rh");
//    depends("t_lapse_rate");
 //   depends("Td_lapse_rate");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}
Walcek_cloud::~Walcek_cloud()
{

};
void Walcek_cloud::run(mesh_elem& face)
{
    //Kunkel RH lapse rates
    // 1/km
    double lapse_rates[] =
            {-0.09,
             0.0,
             0.09,
             0.11,
             0.11,
             0.12,
             0.14,
             0.15,
             0.11,
             0.07,
             -0.02,
             -0.07
            };


    double press_ratio = 0.7;


    double Ta = face->face_data("t");
    double Rh = face->face_data("rh");

//    double Td = mio::Atmosphere::RhtoDewPoint(Rh/100.0,Ta+273.15,false);
//
//
//    double z = face->get_z();
//    double dz = z - 3000.0;// assume 700mb is at 3000m
//
//    double Td_lapse_rate = face->face_data("Td_lapse_rate");
//    double T_lapse_rate = face->face_data("t_lapse_rate");
//
//    double Td_700 = Td + Td_lapse_rate * dz;
//    double Tair_700 = Ta + T_lapse_rate * dz;
//
//    double rh_700 = mio::Atmosphere::DewPointtoRh(Td_700,Tair_700+273.15 ,false);

    double lapse = lapse_rates[global_param->month() - 1] / 1000.0; // -> 1/m
    double rh_700 = Rh * exp(lapse*(3000.0-face->get_z()));

    rh_700 /= 100.0;//factional


    //bound RH
    rh_700 = std::min(1.0,rh_700);
    rh_700 = std::max(0.0,rh_700);

    double dx = 80.0;
    double f_max = 78.0 + dx/15.5; //eqn (2)
    double f_100 = f_max * (press_ratio - 0.1) / 0.6 / 100.0; // eqn (3)
    double one_minus_RHe = 0.196 + (0.76-dx/2834.0) * (1.0 - press_ratio); // eqn (5)


    double cloud_frac = f_100 * exp((rh_700 - 1.0)/one_minus_RHe);
    cloud_frac = std::min(cloud_frac,1.0);


    face->set_face_data("cloud_frac",cloud_frac);

}
