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

//

/**
 * \ingroup modules iswr
 * @{
 * \class fast_shadow
 *
 * Calculates horizon shadows using an adaptation of Dozier and Frew 1990 to an unstructured mesh. Modifies the
 * direct beam iswr ( = 0 W/m^2). Does not impact the diffuse beam
 *
 * **Requires:**
 * - Solar elevation "solar_el" [degrees]
 * - Solar azimuth "solar_az" [degrees]
 *
 * **Provides:**
 * - Binary no shadow/shadow "shadow" [0 or 1=shadow]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "steps": 10,
 *       "max_distance": 1000
 *    }
 *
 *
 * .. confval:: steps
 *
 *    :type: int
 *    :default: 10
 *
 *    Number of steps along the search vector to check for a higher point
 *
 * .. confval:: max_distance
 *
 *    :type: double
 *    :default: 1000 m
 *
 *    Maximum search distance to look for a higher point
 *
 * \endrst
 *
 * **References:**
 * - Dozier, J., & Frew, J. (1990). Rapid calculation of terrain parameters for radiation modeling from digital
 * elevation data. IEEE Transactions on Geoscience and Remote, 28(5), 963–969.
 * - Marsh, C., Pomeroy, J., Wheater, H. (2020). The Canadian Hydrological Model (CHM) v1.0: a multi-scale, multi-extent,
 * variable-complexity hydrological model – design and overview Geoscientific Model Development  13(1), 225-247.
 * https://dx.doi.org/10.5194/gmd-13-225-2020

 * @}
 */
class fast_shadow : public module_base
{
REGISTER_MODULE_HPP(fast_shadow);
public:
    fast_shadow(config_file cfg);

    ~fast_shadow();

    virtual void run(mesh_elem& face);

//number of steps along the search vector to check for a higher point
    int steps;
    //max distance to search
    double max_distance;

    //size of the step to take
    double size_of_step;

};
