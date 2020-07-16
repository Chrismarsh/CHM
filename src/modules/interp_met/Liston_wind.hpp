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
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <cmath>


/**
 * \ingroup modules wind met
 * @{
 * \class Liston_wind
 *
 * Calculates windspeeds using a terrain curvature following Liston and Elder 2006. Direction convention is
 * (North = 0, clockwise from North)
 *
 * **Depends from met:**
 * - Wind at reference height "U_R" [\f$m \cdot s^{-1}\f$]
 * - Direction at reference height 'vw_dir' [degrees] (North = 0, clockwise from North)
 *
 * **Provides:**
 * - Wind  at reference height "U_R"
 * - Interpolated wind field at reference height prior to downscaling "U_R_orig"  [\f$m \cdot s^{-1}\f$]
 * - Wind direction "vw_dir" [degrees]
 * - Original wind direction "vw_dir_orig" [degrees]
 * - Amount the wind vector direction has been changed "vw_dir_divergence" [degrees]
 * - Zonal U at reference height  "zonal_u"  [ \f$ m \cdot s^{-1}\f$ ]
 * - Zonal V at reference height "zonal_v"  [ \f$ m \cdot s^{-1}\f$ ]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "distance": 300,
 *       "Ww_coeff: 1,
 *       "ys": 0.5,
 *       "yc": 0.5
 *    }
 *
 * .. confval:: distance
 *
 *    :type: double
 *    :default: 300
 *
 *    Distance to "look" to compute the terrain curvature.
 *
 * .. confval:: Ww_coeff
 *
 *    :type: int
 *    :default: 1
 *
 *    Leading coefficient in equation 16. Used for compatibility with calibration approach such as done by Pohl.
 *
 * .. confval:: ys
 *
 *    :type: double
 *    :default: 0.5
 *
 *    Slope weight. Valid range [0,1]. The value of 0.5 gives equal weight to slope and curvature
 *
 * .. confval:: yc
 *
 *    :type: double
 *    :default: 0.5
 *
 *    Curvature weight. Valid range [0,1]. The value of 0.5 gives equal weight to slope and curvature
 *
 * \endrst
 *
 * **References:**
 *
 * Liston, G. E., & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet).
 * Journal of Hydrometeorology, 7(2), 217–234. http://doi.org/10.1175/JHM486.1
 * @}
*/
class Liston_wind : public module_base
{
REGISTER_MODULE_HPP(Liston_wind);
public:
    Liston_wind(config_file cfg);
    ~Liston_wind();
    virtual void run(mesh& domain);
    virtual void init(mesh& domain);
    double ys;
    double yc;
    class lwinddata : public face_info
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
    double Ww_coeff;
};
