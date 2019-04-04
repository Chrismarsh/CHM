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

#include "Gray_inf.hpp"
REGISTER_MODULE_CPP(Gray_inf);

Gray_inf::Gray_inf(config_file cfg)
        : module_base("Gray_inf", parallel::data, cfg)
{

    depends("swe");
    depends("snowmelt_int");

    provides("inf");
    provides("total_inf");
    provides("total_excess");
    provides("runoff");
    provides("soil_storage");
    provides("potential_inf");
    provides("opportunity_time");
    provides("available_storage");

}

Gray_inf::~Gray_inf()
{

}

void Gray_inf::init(mesh& domain)
{
    ompException oe;
    //store all of snobals global variables from this timestep to be used as ICs for the next timestep
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
      oe.Run([&]
	     {
	       auto face = domain->face(i);
	       auto* d = face->make_module_data<Gray_inf::data>(ID);

	       d->soil_depth = 400;
	       d->porosity = .4;
	       d->max_storage =  d->soil_depth * d->porosity;
	       //        d->storage =  d->max_storage  - (1 - d->max_storage * face->parameter("sm"_s)/100.);
	       d->storage =  d->max_storage * face->parameter("sm"_s)/100.;

	       d->last_ts_potential_inf = 0;
	       d->opportunity_time=0.;
	       d->total_inf = 0.;
	       d->total_excess = 0.;
	     });
    }
    oe.Rethrow();
}
void Gray_inf::run(mesh_elem &face)
{
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }


    auto* d = face->get_module_data<Gray_inf::data>(ID);

    auto id = face->cell_local_id;
    double C = 2.;
    double S0 = 1;
    double SI = face->get_initial_condition("sm")/100.;

    double TI = 272.;

    double runoff = 0.;
    double inf = 0.;

    double snowmelt = (*face)["snowmelt_int"_s];

    double potential_inf = d->last_ts_potential_inf;
    double avail_storage = (d->max_storage - d->storage);

    if(snowmelt > 0)
    {
        d->opportunity_time += global_param->dt() / 3600.;
        double t0 = d->opportunity_time;
        potential_inf = C * pow(S0,2.92) * pow((1. - SI),1.64) * pow((273.15 - TI) / 273.15, -0.45) * pow(t0,0.44);


        //cap the total infiltration to be no more than our available storage
        if (potential_inf > avail_storage )
        {
            potential_inf = avail_storage;
        }

        d->last_ts_potential_inf = potential_inf;

        if( d->total_inf + snowmelt > potential_inf)
        {
            runoff = (d->total_inf + snowmelt) - potential_inf;
        }


        //infiltrate everything else
        inf    = snowmelt - runoff;


        d->storage += inf;

        if (d->storage > d->max_storage)
        {
            d->storage = d->max_storage;
        }


        d->total_inf += inf;
        d->total_excess += runoff;

    }


    if( !is_nan(face->get_initial_condition("sm")))
    {
        (*face)["total_excess"_s]=d->total_excess;
        (*face)["total_inf"_s]=d->total_inf;

        (*face)["runoff"_s]=runoff;
        (*face)["inf"_s]=inf;
        (*face)["potential_inf"_s]=potential_inf;
        (*face)["soil_storage"_s]= d->storage;
        (*face)["opportunity_time"_s]=d->opportunity_time;
        (*face)["available_storage"_s]=avail_storage;
    }
    else
    {
        set_all_nan_on_skip(face);
    }

}
