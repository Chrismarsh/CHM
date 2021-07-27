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


#include "fetchr.hpp"
REGISTER_MODULE_CPP(fetchr);

fetchr::fetchr(config_file cfg)
        : module_base("fetchr", parallel::data, cfg)
{
    depends("vw_dir");

    provides("fetch");

    //number of steps along the search vector to check for a higher point
    steps = cfg.get("steps",10);
    //max distance to search
    max_distance = cfg.get("max_distance",1000.0);

    //size of the step to take
    size_of_step = max_distance / steps;

    I = cfg.get("I",0.06);

    incl_veg = cfg.get("incl_veg",true);

    h_IBL = 5;

}

fetchr::~fetchr()
{


}

void fetchr::run(mesh_elem& face)
{


    (*face)["fetch"_s]= max_distance;
    return;
    //direction it is from, need upwind fetch
    double wind_dir = (*face)["vw_dir"_s] ;

    //if we are using vegetation and the current face is covered in veg, set the fetch to 0
    if(incl_veg && face->has_vegetation())
    {

        double me_Z_CanTop = face->veg_attribute("CanopyHeight");
        if(me_Z_CanTop > 1) // 1m might be too high?
        {
            (*face)["fetch"_s]= 0;
            return;
        }

    }

    // search along wind_dir azimuth in j step increments
    for (int j = 1; j <= steps; ++j)
    {
        double distance = j * size_of_step;

        auto f = face->find_closest_face(wind_dir, distance);

        double Z_CanTop = 0;
        if (incl_veg && f->has_vegetation())
        {

            Z_CanTop = f->veg_attribute("CanopyHeight");
        }

        //include canopy height if available
        double Z_test =  f->center().z()+Z_CanTop;

        //equation 1, pg 771, Lapen and Martz 1993
        double Z_core = face->center().z() + distance*I;

        double z0_1 = 0.12*Z_CanTop;
        double z0_2 = 0.001;
        double n=1.0/0.8;
        double h = 5; //IBL height, m

        //Kaimal, J., Finnigan, J., 1994. Atmospheric Boundary Layer Flows: Their structure and measurement. Oxford University Press, Toronto.
        //eq 4.2,4.3
        double x_sss  = pow(((33.33333333*h-25.*z0_2)/(log(z0_1/z0_2)*z0_2)),n)*z0_2;

        //reset the fetch if there is elevation, but also if we run into vegetation and haven't restablished the IBL (assumed @ 5m)
        if(Z_test >= Z_core ||
                (incl_veg && distance < x_sss) )
        {
            (*face)["fetch"_s]= distance;
            break;
        }
    }




}
