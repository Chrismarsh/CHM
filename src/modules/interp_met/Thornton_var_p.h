/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "module_base.hpp"
#include <gsl/gsl_fit.h>
#include <vector>
#include <gsl/gsl_combination.h>

/**
* \addtogroup modules
* @{
* \class Precip
* \brief Calculates precipitation
*
* Spatially distributes liquid water precipitation using the Thornton, et al. 1997 method.
* Unlike Thornton_p, this calculates scaling rates based on observed data on a per-timestep basis
*
*
* Depends:
* - Precip from met file "p" [mm]
*
* Provides:
* - Precip "p" [mm]
*
* Reference:
* - Thornton, P. E., Running, S. W., & White, M. A. (1997). Generating surfaces of daily meteorological variables over large regions of complex terrain. Journal of Hydrology, 190(3-4), 214–251. http://doi.org/10.1016/S0022-1694(96)03128-9
* */
class Thornton_var_p : public module_base
{
public:
    Thornton_var_p(config_file cfg);
    ~Thornton_var_p();
    void run(mesh_elem& face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};

