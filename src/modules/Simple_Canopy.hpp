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

#pragma once

#include <boost/shared_ptr.hpp>
#include "logger.hpp"
#include "module_base.hpp"
#include <meteoio/MeteoIO.h>
#include <physics/Atmosphere.h>
#include <physics/PhysConst.h>
#include <physics/Vegetation.h>
#include <string>


/*
 * @Brief Solves energy and mass equations for canopy states using parameterizations based on LAI and canopy closure.
 *
 */

class Simple_Canopy : public module_base
{
REGISTER_MODULE_HPP(Simple_Canopy);
public:
    Simple_Canopy(config_file cfg);

    ~Simple_Canopy();

    virtual void run(mesh_elem &elem);

    virtual void init(mesh& domain);

    double delta(double ta);

    double lambda(double ta);

    double gamma(double air_pressure, double ta);

    double Qs(double air_pressure, double ta);

    struct data : public face_info {

        // parameters
        double LAI;
        double CanopyHeight;
        int canopyType;

        // model states
        double rain_load;
        double Snow_load;

        // Output diagnostic variables
        double cum_net_snow;
        double cum_net_rain;
        double cum_Subl_Cpy;
        double cum_intcp_evap;
        double cum_SUnload_H2O;
    };



};
