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
 * \class iswr
 *
 * This aggregates the _no_slope direct and diffuse flat plane estimates provided by other modules and calculates a unified
 * variant that includes (if required) slope effects and applies shadow masks from shadowing modules.
 *
 * \rst
 * .. warning::
 *    If using ground-based observations, please look at the ``already_cosine_corrected`` configuration option!
 *
 * \endrst
 *
 * **Depends:**
 * - Shortwave direct on a flat plane "iswr_direct_no_slope" [ \f$ W \cdot m^{-2}\f$]
 * - Shortwave diffuse on a flat plane "iswr_diffuse_no_slope" [ \f$ W \cdot m^{-2}\f$]
 *
 * Provides:
 * - Shortwave all beam "iswr" [ \f$ W \cdot m^{-2}\f$]
 * - Shortwave direct "iswr_direct" [ \f$ W \cdot m^{-2}\f$]
 * - Shortwave diffuse "iswr_diffuse" [ \f$ W \cdot m^{-2}\f$]
 *
 * **Configuration:**
 *
 * \rst
 * .. code:: json
 *
 *    {
 *       "no_slope": false,
 *       "already_cosine_corrected": false
 *    }
 *
 *
 * .. confval:: no_slope
 *
 *    :type: boolean
 *    :default: false
 *
 *    Disables slope correction.
 *
 * .. confval:: already_cosine_corrected
 *
 *    :type: boolean
 *    :default: false
 *
 *    "Most observations implicitly have a cosine-correction built in by virtu of the flat-plane observation.
 *     When using observed, ground based observations, you probably want to set this to true, so that a cosine correction
 *     is not applied twice.
 *
 *     If the forcing input is out of a NWP, then this should be ``false``.
 *
 * \endrst
 */
class iswr : public module_base
{
REGISTER_MODULE_HPP(iswr);
    public:
        iswr(config_file cfg);
        ~iswr();
        virtual void run(mesh_elem& face);


        bool assume_no_slope;

        // if we are using obs, then our obs implicitily have a cosine correction.
        // This needs to be undo prior to the correction for slope
        bool already_cosine_corrected;

};

