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
 * \ingroup modules iswr nwp
 * @{
 * \class iswr_from_nwp
 * Spatially interpolates total and diffuse shortwave radiation. Computes direct beam as a residual.
 * Best used from NWP system.
 *
 * **Requires from met:**
 * - Total shortwave radiation met file - "Qsi" [ \f$ W \cdot m^2 \f$]
 * - Total diffuse shortwave radiation met file - "Qsi_diff" [ \f$ W \cdot m^2 \f$]
 *
 * **Provides:**
 * - Shortwave all beam - "iswr" [ \f$ W \cdot m^2 \f$]
 * - Shortwave direct - "iswr_direct" [ \f$ W \cdot m^2 \f$]
 * - Shortwave diffuse - "iswr_diffuse" [ \f$ W \cdot m^2 \f$]
 *
 * **Configuration keys:**
 * - None
 * @}
 */
class iswr_from_nwp : public module_base
{
REGISTER_MODULE_HPP(iswr_from_nwp);
public:
    iswr_from_nwp(config_file cfg);
    ~iswr_from_nwp();
    void run(mesh_elem &face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};
