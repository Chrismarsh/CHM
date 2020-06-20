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


#include "Burridge_iswr.hpp"
REGISTER_MODULE_CPP(Burridge_iswr);

Burridge_iswr::Burridge_iswr(config_file cfg)
        : module_base("Burridge_iswr", parallel::data, cfg)
{

    depends("cloud_frac");
    depends("solar_el"); //degrees

    provides("iswr_diffuse_no_slope"); //Shortwave diffuse beam without slope correction
    provides("iswr_direct_no_slope"); //Shortwave direct beam without slope correction

    provides("atm_trans");
}

Burridge_iswr::~Burridge_iswr()
{

}

void Burridge_iswr::run(mesh_elem &face)
{
    double solar_el = (*face)["solar_el"_s];
    double cosZ = cos( (90.0-solar_el) *mio::Cst::to_rad);

//    double aspect_south0 = face->aspect() * mio::Cst::to_deg;
//    if (aspect_south0 >= 180.0)
//        aspect_south0 -=  180.0;
//    else
//        aspect_south0 += 180.0;
//    aspect_south0 *= mio::Cst::to_rad;
//
//    double slope = face->slope();
//    double sun_az = global_param->solar_az(); //* mio::Cst::to_rad;
//    if (sun_az >= 180.0)
//        sun_az -=  180.0;
//    else
//        sun_az += 180.0;
//    sun_az *= mio::Cst::to_rad;
//
//
//    double sinZ = sqrt(1.0 - cosZ*cosZ);
//    double cosi = cos(slope) * cosZ +
//              sin(slope) * sinZ *
//              cos(sun_az - aspect_south0);
//
//
//    if (cosi < 0.0)
//        cosi = 0.0;
//    if(cosZ <= 0.0)
//        cosZ=0.0;

    double S = 1375.0;
    double cf = (*face)["cloud_frac"_s];

    double dir = S  * (0.6+0.2*cosZ)*(1.0-cf);
    double diff = S * (0.3+0.1*cosZ)*(cf);

 //   dir = dir * cosi;
    diff = diff*cosZ;


    if (diff <0)
        diff = 0.0;
    if(dir <0)
        dir = 0.0;

    (*face)["iswr_diffuse_no_slope"_s]=diff;
    (*face)["iswr_direct_no_slope"_s]=dir;

    //constrain to be [0,1]
    double tau = (dir+diff) / 1375.;
    if(tau < 0)
        tau = 0;
    if(tau > 1)
        tau = 1;

    (*face)["atm_trans"_s] = tau;
}
