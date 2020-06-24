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

/**
 * \ingroup modules snow
 * @{
 * \class Richard_albedo
 *
 * Slow albedo changes for cold snow are neglected, but an exponential decay to an asymptotic minimum of 0.5
 * with an adjustable time constant applied to melting snow.
 *
 * Somewhat mis-named as although this is detailed in Essery and Etchevers (2004), the original source is from CLASS
 * detailed in Verseghy, et al (1991).
 *
 * **Depends:**
 * - Snow Water Equivalent "swe" [mm]
 * - Snow surface temperature "T_s_0" [\f$ K \f$ ]
 * - Precipitation, snow phase "p_snow" [\f$ mm \cdot dt^{-1} \f$ ]
 *
 * **Provides:**
 * - Snow albedo "snow_albedo" [0-1]
 * - True if a melting snow albedo is being used "melting_albedo" [0 or 1]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "albedo_min": 0.5,
 *       "albedo_max": 0.84,
 *       "a1": 1.08e7,
 *       "a2": 7.2e5,
 *       "min_swe_refresh": 1.0,
 *       "init_albedo_snow": 0.85,
 *       "init_albedo_bare": 0.17
 *    }
 * .. confval:: albedo_min
 *
 *    :type: double
 *    :default: 0.5
 *
 *    Minimum snow albedo
 *
 * .. confval:: albedo_max
 *
 *    :type: double
 *    :default: 0.84
 *
 *    Maximum snow albedo
 *
 * .. confval:: a1
 *
 *    :type: double
 *    :default: 1.08e7 s
 *
 *    Cold snow decay constant
 *
 * .. confval:: a2
 *
 *    :type: double
 *    :default: 7.2e5 s
 *
 *    Melting snow decay constant
 *
 * .. confval: min_swe_refresh
 *
 *    :type: double
 *    :default: 1 mm
 *
 *    Minimum precipitation required to fully refresh the snow albedo to max
 *
 * .. confval:: init_albedo_snow
 *
 *    :type: double
 *    :default: 0.85
 *
 *    Initial fresh snow albedo
 *
 * .. confval:: init_albedo_bare
 *
 *    :type: double
 *    :default: 0.17
 *
 *    Bare ground albedo
 * \endrst
 *
 * **References:**
 * - Equation 4 and 5
 * - Essery, R., and P. Etchevers (2004), Parameter sensitivity in simulations of snowmelt, J. Geophys. Res., 109(D20111), 1–15, doi:10.1029/2004JD005036.
 * - Verseghy, D. L.: Class – A Canadian land surface scheme for GCMS, I. Soil model, Int. J. Climatol., 11, 111–133, https://doi.org/10.1002/joc.3370110202, 1991
 * @}
 */
class Richard_albedo : public module_base
{
REGISTER_MODULE_HPP(Richard_albedo);
public:
    struct data : public face_info
    {
        double albedo;

    };

    Richard_albedo(config_file cfg);
    ~Richard_albedo();
    void run(mesh_elem& face);
    void init(mesh& domain);
    void checkpoint(mesh& domain,  netcdf& chkpt);
    void load_checkpoint(mesh& domain,  netcdf& chkpt);

    double amin;
    double amax;
    double a1;
    double a2;
    double albedo_snow;
    double albedo_bare;
    double min_swe_refresh;
};
