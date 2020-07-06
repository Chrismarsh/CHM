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


#include "PenmanMonteith_evaporation.hpp"
REGISTER_MODULE_CPP(PenmanMonteith_evaporation);

PenmanMonteith_evaporation::PenmanMonteith_evaporation(config_file cfg)
        :module_base("PenmanMonteith_evaporation", parallel::data, cfg)
{

    depends("iswr");
    depends("ilwr");
    depends("rh");
    depends("t");
    depends("U_2m_above_srf");

    provides("ET");


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

 void PenmanMonteith_evaporation::run(mesh_elem& face)
{

    double albedo = 0.23; //grass and crops
    double qsi = (*face)["iswr"_s];
    double Lin = (*face)["ilwr"_s];

    double rh = (*face)["rh"_s] / 100.;
    double t = (*face)["t"_s];
    double es = Atmosphere::saturatedVapourPressure(t);
    double ea = rh * es / 1000.; // kpa


    double T = (*face)["t"_s];
    double u = (*face)["U_2m_above_srf"_s];


    double grass_emissivity = 0.9;


    double Qn = (1-albedo)*qsi;
    double sigma = 5.67*pow(10.0,-8.0); //boltzman
    double Lout = sigma * grass_emissivity * pow(T+273,4.0); //assume ground temp = air temp (lol)

    double Rn = Qn + (Lin-Lout);

    double G = 0.1*Rn;

    double delta = ( 4098.0*(0.6108*exp( (17.27*T) / (T+237.3))))/pow(T+237.3,2.0);

    double psy_const = 0.066; //kpa / K

    double latent_heat = 2501.0-2.361*T; //kJ/kg

    double cp = 1.005; //kJ/kg

    double rho = 1.2; //density dry air, take it as const for now.

    double h = 0.01; //veg height

    double z0 = h/7.6; //maybe fix this?

    double kappa = 0.41;

    double ra = pow(log( (10.0-0.67*h)/z0),2.0)/(pow(kappa,2.0)*u); //10cm veg

    double rc = 62.0; //s/m  unstressed

    double E = (delta*(Qn-G)/latent_heat + (rho*cp*(es-ea)/ra))/(delta + psy_const * (1+rc/ra));

    (*face)["ET"_s]= E;


}

PenmanMonteith_evaporation::~PenmanMonteith_evaporation()
{

}
