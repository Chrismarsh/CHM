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
#include "math/coordinates.hpp"
#include <physics/Atmosphere.h>
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <cmath>


#include <viennacl/linalg/gmres.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>

/**
 * \ingroup modules met wind
 * @{
 * \class MS_wind
 * Calculates windspeeds using the Mason Sykes wind speed from Essery, et al (1999), using 8 lookup map from DBSM.
 * Optionally uses Ryan, et al. for wind direction correction.
 *
 * **Depends from met:**
 * - Wind at reference height "U_R" [ \f$ m \cdot s^{-1}\f$ ]
 * - Direction at reference height "vw_dir" [degrees]
 *
 * **Provides:**
 * - Wind speed at reference height "U_R" [ \f$ m \cdot s^{-1}\f$ ]
 * - Wind direction 'vw_dir' at reference height [degrees]
 * - Speedup magnitude "W_speedup" [-]
 * - Zonal u speed component at 2m "2m_zonal_u" [ \f$ m \cdot s^{-1}\f$ ]
 * - Zonal v speed component at 2m "2m_zonal_v" [ \f$ m \cdot s^{-1}\f$ ]
 * - Original interpolated wind direction "vw_dir_orig" [degrees]
 *
 * **Parameters:**
 *
 * Requires speedup, u, and v parameters named "MS%i_U" and "MS%i_V" and "MS%i" for each of the 8 directions.
 * These need to be generated using DBSM in ``tools/MSwind``.
 *
 * **Configuration:**
 * \rst
 * .. code::
 * 
 *    {
 *       "speedup_height": 2.0,
 *       "use_ryan_dir", false
 *    }
 *
 * .. confval:: speedup_height
 *
 *    :type: double
 *    :default: 2.0
 *
 *    The height at which the MS Wind tool was run for. The default is 2 m and shouldn't be changed.
 *
 * .. confval:: use_ryan_dir
 *
 *    :type: boolean
 *    :default: false
 *
 *    Instead of using the _u and _v components to compute direction perturbation, use the algorithm of Ryan as per Liston and Elder (2006)
 * \endrst
 *
 * **Reference:**
 *
 * Essery, R., Li, L., Pomeroy, J. (1999). A distributed model of blowing snow over complex terrain
 * Hydrological Processes  13(), 2423-2438.
 *
 * @}
 */
class MS_wind : public module_base
{
REGISTER_MODULE_HPP(MS_wind);
public:
    MS_wind(config_file cfg);
    ~MS_wind();
    virtual void run(mesh& domain);
    virtual void init(mesh& domain);
    double ys;
    double yc;
    class data : public face_info
    {
    public:
        double curvature;
        interpolation interp;
        double corrected_theta;
        double W;
        double temp_u;
        interpolation interp_smoothing;
    };
    double distance;
    bool use_ryan_dir;
    double speedup_height; // height at which the speedup is for
};

