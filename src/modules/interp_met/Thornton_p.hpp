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

#include <constants/Atmosphere.h>

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
* \addtogroup modules
* @{
* \class Precip
* \brief Calculates precipitation
*
* Spatially distributes liquid water precipitation using Thornton, et al. 1997. Monthly scaling factors are used, as
 * summarized in Liston, 2006.
*
* Depends:
* - Precip from met file "p" [mm]
*
* Provides:
* - Precip "p" [mm]
* - Precip "p_no_slope" [mm]
* 
 * Reference:
 * - Thornton, P. E., Running, S. W., & White, M. A. (1997). Generating surfaces of daily meteorological variables over large regions of complex terrain. Journal of Hydrology, 190(3-4), 214–251. http://doi.org/10.1016/S0022-1694(96)03128-9
 * - Liston, G. E., & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). Journal of Hydrometeorology, 7(2), 217–234. http://doi.org/10.1175/JHM486.1
*/
class Thornton_p : public module_base
{
REGISTER_MODULE_HPP(Thornton_p);
public:
    Thornton_p(config_file cfg);
    ~Thornton_p();
    virtual void run(mesh_elem& face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };

    // Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
    bool apply_cosine_correction;

};

/**
@}
*/
