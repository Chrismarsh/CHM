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
 * \ingroup modules met precip
 * @{
 * \class Thornton_p
 * Spatially distributes liquid water precipitation using Thornton, et al. 1997 via monthly lapse rates.
 *
 * **Depends from met:**
 * - Precipitation  "p" [\f$mm \cdot dt^{-1}\f$]
 *
 * Provides:
 * - Lapsed precipitation "p" [\f$mm \cdot dt^{-1}\f$]
 * - Precipitation corrected for triangle slope. If ``"apply_cosine_correction": false``, then no change. "p_no_slope" [\f$mm \cdot dt^{-1}\f$]
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
 * \endrst
 *
 * **References:**
 * - Thornton, P. E., Running, S. W., & White, M. A. (1997). Generating surfaces of daily meteorological variables over large regions of complex terrain. Journal of Hydrology, 190(3-4), 214–251. http://doi.org/10.1016/S0022-1694(96)03128-9
 *
 * @}
 */
class Thornton_p : public module_base
{
REGISTER_MODULE_HPP(Thornton_p);
public:
    Thornton_p(config_file cfg);
    ~Thornton_p();
    virtual void run(mesh_elem& face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };

    // Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
    bool apply_cosine_correction;

};
