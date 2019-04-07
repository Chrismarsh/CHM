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

#include <cstdlib>
#include <string>
#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

#include <boost/math/tools/roots.hpp>

//#include <meteoio/MeteoIO.h>
/**
* \addtogroup modules
* @{
* \class Harder_precip_phase
* \brief Calculates precip phase based
*
* Calculates precipitation phase via falling hydrometeor energy balance
*
* Depends:
* - Air temperature "t" [C]
* - Relative Humidity 'rh' [C]
* - Precip "p" [mm]
*
* Provides:
* - Snow precip p_snow [mm]
* - Liquid precip p_rain [mm]
* - Fractional rain frac_precip_rain [-]
* - Fractional snow frac_precip_snow [-]
* - Cumulated Snow precip p_snow [mm]
* - Cumulated Liquid precip p_rain [mm]
*/
class Harder_precip_phase : public module_base
{
REGISTER_MODULE_HPP(Harder_precip_phase);
public:
    Harder_precip_phase(config_file cfg);
    ~Harder_precip_phase();
    virtual void run(mesh_elem& face);
    void init(mesh& domain);
    double b;
    double c;

    class data : public face_info
    {
    public:
        double hours_since_snowfall;
        double acc_rain;
        double acc_snow;

    };


};

/**
@}
*/
