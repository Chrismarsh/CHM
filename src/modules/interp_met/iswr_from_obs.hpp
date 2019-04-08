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
#include <math.h>
#include <algorithm>
#include <meteoio/MeteoIO.h>


/**
* \addtogroup modules
* @{
* \class iswr_from_obs
* \brief Spatially interpolates observed shortwave measurements
*
* Spatially interpolates observed shortwave measurements and splits it into direct and diffuse beams. This split
 * is via Iqbal
*
* Depends:
* -  Shortwave radiation met file "Qsi" [W/m^2]
*
* Provides:
* - Shortwave all beam "iswr" [W/m^2]
* - Shortwave direct "iswr_direct" [W/m^2]
* - Shortwave diffuse "iswr_diffuse" [W/m^2]
* - Atmospheric transmittance, [0,1] "atm_trans" [-]
*/
class iswr_from_obs : public module_base
{
REGISTER_MODULE_HPP(iswr_from_obs);
public:
    iswr_from_obs(config_file cfg);
    ~iswr_from_obs();
    void run(mesh_elem &face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};
