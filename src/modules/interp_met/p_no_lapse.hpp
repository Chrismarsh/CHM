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

#include <physics/Atmosphere.h>

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>
#include <gsl/gsl_fit.h>
/**
 * \ingroup modules precip met
 * @{
 * \class p_no_lapse
 * Spatially distributes liquid water precipitation with no lapsing
 *
 * **Depends from met**:
 * - Precipitation "p" [\f$mm \cdot dt^{-1}\f$]
 *
 * **Provides:**
 * -  Precipitation "p" [\f$mm \cdot dt^{-1}\f$]
 * - Precipitation corrected for triangle slope. If ``"apply_cosine_correction": false``, then no change.  "p_no_slope" [\f$mm \cdot dt^{-1}\f$]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "apply_cosine_correction": false
 *    }
 *
 *
 * .. confval:: apply_cosine_correction
 *
 *    :type: boolean
 *    :default: false
 *
 *    Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
 *
 * \endrst
 * @}
*/
class p_no_lapse : public module_base
{
REGISTER_MODULE_HPP(p_no_lapse);
public:
    p_no_lapse(config_file cfg);
    ~p_no_lapse();
    virtual void run(mesh_elem& face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };

    // Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
    bool apply_cosine_correction;

};

