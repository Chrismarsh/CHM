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
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

/**
 * \ingroup modules iswr
 * @{
 * \class slope_iswr
 *
 * Calculates incoming direct-beam shortwave solar radiation to slope
 *
 * **Depends:**
 * - Incoming solar shortwave radiation, all beam "iswr"  [\f$ W \cdot m^{-2} \f$ ]
 * - Incoming solar shortwave radiation, direct beam "iswr_direct"  [\f$ W \cdot m^{-2} \f$ ]
 * - Incoming solar shortwave radiation, diffuse beam "iswr_diffuse"  [\f$ W \cdot m^{-2} \f$ ]
 * - Solar elevation "solar_el" [degrees]
 * - Solar azimuth "solar_az" [degrees]
 *
 * **Optional:**
 * - Shadow map "shadow" [0,1]
 *
 * **Provides:**
 * - Slope corrected shortwave all beam "iswr" [W/m^2]
 * - Slope corrected shortwave direct "iswr_direct" [W/m^2]
 * - Solar angle "solar_angle" [degrees]
 *
 * \rst
 * .. code:: json
 *
 *    {
 *       "no_slope": false
 *    }
 *
 * .. confval:: no_slope
 *
 *    :default: false
 *
 *    Disables slope correction calcuation when ``true``
 *
 * \endrst
 *
 * @}
 */
class slope_iswr : public module_base
{
REGISTER_MODULE_HPP(slope_iswr);
    public:
        slope_iswr(config_file cfg);
        ~slope_iswr();
        virtual void run(mesh_elem& face);


        bool assume_no_slope;
};

