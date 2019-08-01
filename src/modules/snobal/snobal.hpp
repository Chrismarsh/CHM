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
#include "TPSpline.hpp"
#include <meteoio/MeteoIO.h>

#include "sno.h"
#include "snomacros.h"
class snodata : public face_info
{
public:
    sno data;
    double sum_runoff;
    double sum_melt;
    int dead;
    double delta_avalanche_snowdepth;
    double delta_avalanche_swe;

};
class snobal : public module_base
{
REGISTER_MODULE_HPP(snobal);
public:
    snobal(config_file cfg);

    ~snobal();

    double drift_density; // if we have blowing snow, this is the density of those particles
    double const_T_g; // constant ground temp, degC

    bool use_slope_SWE; // use a slope corrected SWE for compaction eqn

    virtual void run(mesh_elem &face);
    virtual void init(mesh& domain);
    void checkpoint(mesh& domain, netcdf& chkpt);
    void load_checkpoint(mesh& domain, netcdf& chkpt);

};
