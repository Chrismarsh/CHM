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

#include "Harder_precip_phase.hpp"
REGISTER_MODULE_CPP(Harder_precip_phase);

Harder_precip_phase::Harder_precip_phase(config_file cfg)
        :module_base("Harder_precip_phase", parallel::data, cfg)
{
    depends("t");
    depends("rh");
    depends("p");

    provides("Ti");
    provides("frac_precip_rain");
    provides("frac_precip_snow");
    provides("p_snow");
    provides("p_rain");

    provides("acc_snow");
    provides("acc_rain");


    provides("p_snow_hours"); // hours since snowfall

    //     default values:
    //     b=2.630006;
    //     c=0.09336;
    b = cfg.get("const.b",2.630006);
    c = cfg.get("const.c",0.09336);



    LOG_DEBUG << "Successfully instantiated module " << this->ID;


}
Harder_precip_phase::~Harder_precip_phase()
{

}
void Harder_precip_phase::init(mesh& domain)
{

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<data>(ID);
        d.hours_since_snowfall = 0;
        d.acc_rain = 0;
        d.acc_snow = 0;
    }

}
void Harder_precip_phase::run(mesh_elem& face)
{
    double Ta = (*face)["t"_s]+273.15; //K
    double T =  (*face)["t"_s];
    double RH = (*face)["rh"_s];
    double ea = RH/100 * 0.611*exp( (17.3*T) / (237.3+T));

    // (A.6)
    double D = 2.06 * pow(10,-5) * pow(Ta/273.15,1.75);

    // (A.9)
    double lambda_t = 0.000063 * Ta + 0.00673;

    // (A.10) (A.11)
    double L;
    if(T < 0.0)
    {
        L = 1000.0 * (2834.1 - 0.29 *T - 0.004*T*T);
    }
    else
    {
        L = 1000.0 * (2501.0 - (2.361 * T));
    }

    /*
     * The *1000 and /1000 are important unit conversions. Doesn't quite match the harder paper, but Phil assures me it is correct.
     */
    double mw = 0.01801528 * 1000.0; //[kgmol-1]
    double R = 8.31441 /1000.0; // [J mol-1 K-1]

    double rho = (mw * ea) / (R*Ta);

    auto fx = [=](double Ti)
    {
        return boost::math::make_tuple(
                T+D*L*(rho/(1000.0)-.611*mw*exp(17.3*Ti/(237.3+Ti))/(R*(Ti+273.15)*(1000.0)))/lambda_t-Ti,
                D*L*(-0.6110000000e-3*mw*(17.3/(237.3+Ti)-17.3*Ti/pow(237.3+Ti,2))*exp(17.3*Ti/(237.3+Ti))/(R*(Ti+273.15))+0.6110000000e-3*mw*exp(17.3*Ti/(237.3+Ti))/(R*pow(Ti+273.15,2)))/lambda_t-1);
    };

    double guess = T;
    double min = -50;
    double max = 50;
    double digits = 6;

    double Ti = boost::math::tools::newton_raphson_iterate(fx, guess, min, max, digits);

    double frTi = 1.0 / (1.0+b*pow(c,Ti));

    frTi = std::trunc(100.0*frTi) / 100.0; //truncate to 2 decimal positions

    (*face)["Ti"_s]=Ti;
    (*face)["frac_precip_rain"_s]=frTi;
    (*face)["frac_precip_snow"_s]=1.0-frTi;

    double p = (*face)["p"_s];

    (*face)["p_rain"_s]= p * frTi;
    (*face)["p_snow"_s]= p * (1.0-frTi);

    auto& d = face->get_module_data<data>(ID);
    if( p * (1.0-frTi) > 0) // it's snowing
    {
        d.hours_since_snowfall = 0; // reset
    }
    else
    {
        d.hours_since_snowfall  +=  (global_param->dt() / 3600.0) ; // dt(s) -> hr
    }

    (*face)["p_snow_hours"_s]=d.hours_since_snowfall;

    d.acc_rain += p * frTi;
    d.acc_snow += p * (1.0-frTi);

    (*face)["acc_rain"_s]=d.acc_rain;
    (*face)["acc_snow"_s]=d.acc_snow;



}

void Harder_precip_phase::checkpoint(mesh& domain,  netcdf& chkpt)
{

    chkpt.create_variable1D("Harder_precip_phase:hours_since_snowfall", domain->size_faces());
    chkpt.create_variable1D("Harder_precip_phase:acc_rain", domain->size_faces());
    chkpt.create_variable1D("Harder_precip_phase:acc_snow", domain->size_faces());

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        chkpt.put_var1D("Harder_precip_phase:hours_since_snowfall",i,face->get_module_data<data>(ID).hours_since_snowfall);
        chkpt.put_var1D("Harder_precip_phase:acc_rain",i,face->get_module_data<data>(ID).acc_rain);
        chkpt.put_var1D("Harder_precip_phase:acc_snow",i,face->get_module_data<data>(ID).acc_snow);
    }

}
void Harder_precip_phase::load_checkpoint(mesh& domain, netcdf& chkpt)
{

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->get_module_data<data>(ID).hours_since_snowfall = chkpt.get_var1D("Harder_precip_phase:hours_since_snowfall",i);
        face->get_module_data<data>(ID).acc_rain = chkpt.get_var1D("Harder_precip_phase:acc_rain",i);
        face->get_module_data<data>(ID).acc_snow = chkpt.get_var1D("Harder_precip_phase:acc_snow",i);

        (*face)["acc_rain"_s]=face->get_module_data<data>(ID).acc_rain;
        (*face)["acc_snow"_s]=face->get_module_data<data>(ID).acc_snow;
    }


}
