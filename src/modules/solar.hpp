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

#include "module_base.hpp"
#include <ogr_spatialref.h>


/**
 * \ingroup modules iswr
 * @{
 * \class solar
 * Calculates solar position
 *
 * **Depends:**
 * - None
 *
 * **Provides:**
 * - Solar elevation "solar_el" [degrees]
 * - Solar azimuth "solar_az" [degrees]
 *
 * **Provides:**
 * - Sky view view "svf" [-]
 *
 * **Configuration:**
 *
 * \rst
 * .. code:: json
 *
 *    {
 *       "svf":
 *       {
 *          "steps": 10.
 *          "max_distance": 1000.0,
 *          "nsectors": 12,
 *          "compute": true
 *       }
 *    }
 *
 * .. confval:: steps
 *
 *    :type: int
 *    :default: 10
 *
 *    Number of search steps
 *
 * .. confval:: max_distance
 *
 *    :default: 1000.0
 *
 *    Maximum search distance
 *
 * .. confval:: nsectors
 *
 *    :type: int
 *    :default: 12
 *
 *    Number of azimuthal search bins
 *
 * .. confval:: compute
 *
 *    :default: true
 *
 *    Compute the sky view factor
 *
 *
 * \endrst
 *
 * @}
 */
class solar : public module_base
{

REGISTER_MODULE_HPP(solar);
public:

    //if we have a UTM mesh, cache the calculated lat and long
    struct data : public face_info
    {
        double lat;
        double lng;
    };

    solar(config_file cfg);
    ~solar();
    void run(mesh_elem &face);
    void init(mesh& domain);
};
