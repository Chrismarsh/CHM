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


#include "Evapotranspiration_All.hpp"

REGISTER_MODULE_CPP(Evapotranspiration_All);

Evapotranspiration_All::Evapotranspiration_All(config_file cfg)
        :module_base("Evapotranspiration_All", parallel::data, cfg)
{
    // TODO Constructor is not properly editted with all new inputs (see set vars function at the end)
    depends("iswr");
    depends("netall");
    depends("P_atm");
    depends("ea");
    depends("t");
    depends("U_2m_above_srf"); // 
    depends("soil_storage");                      // but is how albedo is used in CHM as of Sept, 2024

    provides("ET");
    provides("stomatal_resistance");

}

void Evapotranspiration_All::init(mesh& domain)
{
    alpha = cfg.get("alpha_PriestleyTaylor",1.26);
    wind_height = cfg.get("wind_measurement_height",2);
    stomatal_resistance_min = cfg.get("stomatal_resistance_min",62);
    Frac_to_ground = cfg.get<double>("Frac_to_ground",1.0); // iswr_subcanopy exists.

    SoilDataObj = std::make_unique<Soil::soils_na>();

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto& d = face->make_module_data<Evapotranspiration_All::data>(ID);
        
        // Consider if an if statement is necessary.

        
        d.soil_depth = face->soil_attribute<double>("soil_depth"_s);
        if (face->has_vegetation())
        {
            d.LAI = face->veg_attribute("LAI");
            d.LAImax = face->veg_attribute("LAImax");
            d.vegetation_height = face->veg_attribute("CanopyHeight");
        }
        else
        {
            d.LAI = 0.0;
            d.LAImax = 0.0;
            d.vegetation_height = 0.0;
        }
        // TODO Consider if we can have a single model object per triangle.
        // Polymorpish would let this work well
        // It wouldn't work if PT model is used for normal soils when they are saturated.
        if (is_water(face))
        {
            init_PriestleyTaylor(d,alpha);
        }
        // TODO add wetlands
        else 
        {
            init_PenmanMonteith(d, face, wind_height, stomatal_resistance_min, Frac_to_ground);
        }
        
    }

}

void Evapotranspiration_All::run(mesh_elem& face)
{
    auto& d = face->get_module_data<Evapotranspiration_All::data>(ID);

    if (is_water(face))
    {
        // Do PriestlyTaylor
        PT_vars my_PT_vars = set_PriestleyTaylor_vars(face);
        model_output output;
        d.MyPriestleyTaylor->CalcEvapT(my_PT_vars,output);
        
        (*face)["ET"_s] = output.ET;
    }
    // TODO wetlands or saturated soils, need input from soils module.
    // else if (If wetlands)
    // {
    //     Also Do PriestlyTaylor
    // }
    else
    {
        // All members of PM_vars are references and must be set at initialization
        // Therefore we need a copy of SVP to reference
        // t is made its own copy to avoid dereferencing face for "t" twice

        double t = (*face)["t"_s];
        double SVP = Atmosphere::saturatedVapourPressure(t+273.15)/1000; // units of kelvin expected 
        double VP = (*face)["ea"_s];
        PM_vars my_PM_vars = set_PenmanMonteith_vars(face,t,SVP,VP);
        PM_output output;
        d.MyPenmanMonteith->CalcEvapT(my_PM_vars,output);

        (*face)["stomatal_resistance"_s] = output.stomatal_resistance;

        (*face)["ET"_s] = output.ET;
    }
    
    // TODO total_ET, as well as PT ET and PM ET as separate. 
}

Evapotranspiration_All::~Evapotranspiration_All()
{

}

void Evapotranspiration_All::init_PriestleyTaylor(Evapotranspiration_All::data& d,double& alpha)
{
    d.MyPriestleyTaylor = std::make_unique<PriestleyTaylor>(alpha,Atmosphere::Cp);
}

void Evapotranspiration_All::init_PenmanMonteith(Evapotranspiration_All::data& d,mesh_elem& face, double& wind_height, double& stomatal_resistance_min, double& Frac_to_ground)
{
    
    const std::string soil_type = face->soil_attribute<std::string>("soil_type"_s,"soils");
    
    const double Cp = Atmosphere::Cp;
    const double kappa = Atmosphere::kappa;
    const double air_entry_tension = SoilDataObj->air_entry_tension(soil_type);
    const double pore_size_dist = SoilDataObj->pore_size_dist(soil_type);
    const double wilt_point = SoilDataObj->wilt_point(soil_type);
    const double porosity = SoilDataObj->porosity(soil_type);
    
    double soil_storage_max = face->soil_attribute<double>("soil_storage_max"_s);
    // Leaf area index is not used if no vegetation, but LAI and LAImax are references in the PenmanMonteith model, therefore values are needed for initialization. It is ok if these values go out of scope as long as there is no vegetation. 

    d.MyPenmanMonteith = std::make_unique<PenmanMonteith>(d.LAI, d.LAImax, d.vegetation_height, wind_height, 
            stomatal_resistance_min, d.soil_depth, Frac_to_ground, Cp, kappa, 
            air_entry_tension, pore_size_dist, wilt_point, porosity);
    
}

PM_vars Evapotranspiration_All::set_PenmanMonteith_vars(mesh_elem& face,double& t, double& saturated_vapour_pressure,double& vapour_pressure)
{
    PM_vars vars((*face)["U_2m_above_srf"_s],(*face)["iswr"_s],(*face)["netall"_s],t,(*face)["soil_storage"_s],vapour_pressure,saturated_vapour_pressure,(*face)["P_atm"_s]); 
    
    return vars;
}

PT_vars Evapotranspiration_All::set_PriestleyTaylor_vars(mesh_elem& face)
{
    // TODO P_atm, is a state variable and should have a copy local to the run function 
    // because it is calculated from the Atmosphere namespace, not done yet
    PT_vars vars((*face)["netall"_s],(*face)["P_atm"_s],(*face)["t"_s]);

    return vars;
}
