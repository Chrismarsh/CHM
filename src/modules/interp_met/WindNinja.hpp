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

#include <Winstral_parameters.hpp>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <cmath>


/**
 * \ingroup modules met wind
 * @{
 * \class WindNinja
 *
 * Calculates wind speed and direction following the downscaling stategy of Barcons et al. (2018). This is via a library
 * of high-resolution wind field generated with the WindNinja wind flow model.
 *
 * **Depends from met:**
 * - Wind at reference height "U_R" [ \f$ m \cdot s^{-1}\f$ ]
 * - Direction at reference height "vw_dir" [degrees]
 *
 * **Provides:**
 * - Wind speed at reference height "U_R" [ \f$ m \cdot s^{-1}\f$ ]
 * - Wind direction 'vw_dir' at reference height [degrees]
 * - Zonal U at reference height  "zonal_u"  [ \f$ m \cdot s^{-1}\f$ ]
 * - Zonal V at reference height "zonal_v"  [ \f$ m \cdot s^{-1}\f$ ]
 *
 * **Diagnostics:**
 * - Downscaled windspeed "Ninja_speed" [ \f$ m \cdot s^{-1}\f$ ]
 * - Interpolated wind field at reference height prior to downscaling "U_R_orig"  [\f$m \cdot s^{-1}\f$]

 * - Interpolated zonal U at level H_Forc "interp_zonal_u"  [ \f$ m \cdot s^{-1}\f$ ]
 * - Interpolated zonal V at level H_Forc "interp_zonal_v"  [ \f$ m \cdot s^{-1}\f$ ]
 * - What lookup map was used "lookup_d"  [-]
 * - Original wind direction "vw_dir_orig" [degrees]
 * - Amount the wind vector direction has been changed "vw_dir_divergence" [degrees]
 *
 * **Configuration:**
 * \rst
 * .. code::
 *
 *    {
 *       "ninja_average": true,
 *       "compute_Sx": true,
 *       "ninja_recirc": false,
 *       "Sx_crit", 30.
 *       "L_avg": 1000.
 *       "H_forc": 40,
 *       "Max_spdup": 3.0,
 *       "Min_spdup": 0.1,
 *    }
 *
 *
 * .. confval:: ninja_average
 *
 *    :type: boolean
 *    :default: true
 *
 *    Linear interpolation between the closest 2 wind fields from the library
 *
 * .. confval:: compute_Sx
 *
 *    :type: boolean
 *    :default: true
 *
 *   Use the Winstral Sx parameterization to idenitify and modify lee-side windfield. This will cause a runtime error (conflict) if
 *   ``Winstral_parameters`` is also a module. This uses an angular window of 30 degrees and a step size of 10 m.
 *
 * .. confval:: ninja_recirc
 *
 *    :type: boolean
 *    :default: false
 *
 *    Enables the leeside slow down via ``compute_Sx``. Requires ``"compute_Sx":true``.
 *
 * .. confval:: Sx_crit
 *
 *    :type: double
 *    :default: 30.
 *
 *    Reduce wind speed on the lee side of mountain crest identified by Sx>Sx_crit
 *
 * ..confval:: L_avg
 *
 *    :type: int
 *    :default: None
 *
 *    The WindMapper tool uses a radius to compute a mean. This is the length over which that average is done. Normally this will be
 *    baked into the parameter name (e.g., Ninja_1000). However, there may be reasons to specifiy it directly.
 *    Cannot be used if the parameters have the _Lavg suffix.
 *    You likely don't need to set this.
 *
 * .. confval:: H_forc
 *
 *    :type: double
 *    :default: 40.0
 *
 *     Reference height for input forcing and WindNinja wind field library
 *
 * .. confval:: Max_spdup
 *
 *    :type: double
 *    :default: 3.0
 *
 *    Limit speed up value to Max_spdup to avoid unrelistic values at crest top
 *
 * .. confval:: Min_spdup
 *
 *    :type: double
 *    :default: 3.0
 *
 *    Limit speed up value to Max_spdup to avoid unrelistic values at crest top
 *
 * \endrst
 *
 * **Parameters:**
 *
 * Requires speedup, u, and v parameters named "Ninja%i_U" and "Ninja%i_V" and "Ninja%i" for each of the _n_ directions.
 * Should be generated with WindMapper. The number of directions will be automatically determined as will the Lavg value.
 * These should be computed with the <a href=https://windmapper.readthedocs.io/en/latest/index.html> Windmapper tool </a>.
 *
 *
 * **References:**
 * - Barcons, J., Avila, M., Folch, A. (2018). A wind field downscaling strategy based on domain segmentation and transfer functions Wind Energy  21(6)https://dx.doi.org/10.1002/we.2169
 *
 * @}
 */
class WindNinja : public module_base
{
REGISTER_MODULE_HPP(WindNinja);
public:
    WindNinja(config_file cfg);
    ~WindNinja();
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
        double W_transf;
    };
    double distance;
    int N_windfield; //  Number of wind fields in the library
    bool ninja_average; // Boolean to activate linear interpolation betweem the closest 2 wind fields from the library
    double H_forc; // Reference height for GEM forcing and WindNinja wind field library
    double Max_spdup;  // Maximal value of crest speedup
    double Min_spdup;  // Minimal value of crest speedup
    int L_avg;   // Size of the square of averaging when computing the speed-up map
                 // Not used by default and set to -1 if not specified in the config file.
    bool ninja_recirc; // Boolean to activate wind speed reduction on the leeside of mountainous terrain

    bool compute_Sx; // uses the Sx module to influence the windspeeds so Sx needs to be computed during the windspeed evaluation, instead of a seperate module
    double Sx_crit;    // Critical values of the Winstral parameter to determine the occurence of flow separation.
    boost::shared_ptr<Winstral_parameters> Sx;
};
