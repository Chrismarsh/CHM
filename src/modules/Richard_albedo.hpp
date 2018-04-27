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

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <meteoio/MeteoIO.h>

/**
 * Eqn 4, 5
 * Essery, R., and P. Etchevers (2004), Parameter sensitivity in simulations of snowmelt, J. Geophys. Res., 109(D20111), 1–15, doi:10.1029/2004JD005036.
 */
class Richard_albedo : public module_base
{
public:
    struct data : public face_info
    {
        double albedo;

    };

    Richard_albedo(config_file cfg);
    ~Richard_albedo();
    void run(mesh_elem& face);
    void init(mesh domain);
    void checkpoint(mesh domain,  netcdf& chkpt);
    void load_checkpoint(mesh domain,  netcdf& chkpt);

    double amin;
    double amax;
    double a1;
    double a2;
    double albedo_snow;
    double albedo_bare;
    double min_swe_refresh;
};


