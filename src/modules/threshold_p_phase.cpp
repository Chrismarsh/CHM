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

#include "threshold_p_phase.hpp"
REGISTER_MODULE_CPP(threshold_p_phase);

threshold_p_phase::threshold_p_phase(config_file cfg)
    : module_base("threshold_p_phase", parallel::data, cfg)
{
    depends("t");
    depends("p");

    t_thresh = cfg.get("threshold_temperature", 2.0) ;

    provides("frac_precip_rain");
    provides("frac_precip_snow");

    provides("p_snow");
    provides("p_rain");

}

threshold_p_phase::~threshold_p_phase()
{

}

void threshold_p_phase::run(mesh_elem &face)
{
    double p = (*face)["p"_s];

    if( (*face)["t"_s] >= t_thresh)
    {
        (*face)["frac_precip_rain"_s]= 1;
        (*face)["frac_precip_snow"_s]= 0;

        (*face)["p_rain"_s]= p;
        (*face)["p_snow"_s]= 0 ;
    } else
    {
        (*face)["frac_precip_rain"_s]= 0;
        (*face)["frac_precip_snow"_s]= 1;

        (*face)["p_rain"_s]= 0;
        (*face)["p_snow"_s]= p ;
    }



}
