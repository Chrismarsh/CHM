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

#include <math.h>
#include <constants/Snow.h>

#pragma once

namespace Atmosphere {
    /********* Atmosphere ************/

    const double Z_U_R = 50.0; // Reference height all input wind speed forcing data are scaled to so that 1) they can be interpolated and 2) all modules can scale
    // down the reference wind speed to the height they require (i.e. through a canopy)

    const double KinVisc = 1.88e-5; // kinematic viscosity of air (Sask. avg. value) (units ????)

    double log_scale_wind(double u, double Z_in, double Z_out, double snowdepthavg, double z0=Snow::Z0_SNOW);

    // Inoue E (1963) On the turbulent structure of air flow within crop canopies. J Meteorol Soc Jpn 41:317–326
    double exp_scale_wind(double u, double Z_in, double Z_out, const double alpha);

   // Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
   // See Fig 1 and 2 of Kienzle (2010, Hydrological Processes)
    double corr_precip_slope(double p, double slope);


}


