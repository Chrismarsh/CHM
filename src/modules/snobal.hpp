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
#include <meteoio/MeteoIO.h>

#include "snobal/sno.h"
#include "snobal/snomacros.h"
class snodata : public face_info
{
public:
    sno data;
    double sum_runoff;
    double sum_subl;
    double sum_pcp_sno;
    double sum_melt;
    int dead;
    double delta_avalanche_snowdepth;
    double delta_avalanche_swe;

};

/**
 * \ingroup modules snow
 * @{
 * \class snobal
 *
 * Snobal is a physically-based snowpack model that approximates the snowpack with two layers. The surface-active layer
 * has a fixed thickness of 0.1 m and is used to estimate surface temperature for outgoing longwave radiation and turbulent heat fluxes.
 * The second lower layer represents the remaining snowpack. For each layer, Snobal simulates the evolution of the snow water equivalent,
 * temperature, density, cold content, and liquid water content. This version of Snobal includes an improved algorithm for
 * snow compaction that accounts for bulk compaction and temperature metamorphism (Hedrick et al., 2018).
 *
 * **Depends:**
 * - Incoming shortwave radiation, all beams "iswr" [ \f$ W \cdot m^{-2} \f$ ]
 * - Incomging longwave radiation "ilwr" [ \f$ W \cdot m^{-2} \f$ ]
 * - Relative humidity "rh" [%]
 * - Air temperature "t" [\f$ {}^\circ C \f$]
 * - Windspeed at 2m "U_2m_above_srf" [\f$ m \cdot s^{-1}    \f$]
 * - Precipitation "p" [\f$ mm \cdot dt^{-1} \f$]
 * - Precipitation snow fraction "frac_precip_snow" [-]
 * - Snow albedo "snow_albedo" [-]
 *
 * **Optional:**
 * - Optionally, depend on the ``_subcanopy`` variants of the above
 * - Blowing snow erosion/deposition mass "drift_mass" [mm]
 * - Ground temperature "T_g" [\f$  {}^\circ C   \f$]
 * - Change in snow mass due to avalanching "delta_avalanche_mass" [mm]
 * - Change in snow depth due to avalanching "delta_avalanche_snowdepth" [mm]
 *
 * **Provides:**
 * - Snow water equivalent "swe" [mm]
 * - Interval snowmelt "snowmelt_int" [mm]
 * - Net allwave radiation "R_n" [ \f$ W \cdot m^{-2} \f$ ]
 * - Sensible heat flux "H" [ \f$ W \cdot m^{-2} \f$ ]
 * - Latent heat flux "E" [ \f$ W \cdot m^{-2} \f$ ]
 * - Ground heat flux "G" [ \f$ W \cdot m^{-2} \f$ ]
 * - Advected heat from precipitation "M" [ \f$ W \cdot m^{-2} \f$ ]
 * - Change in internal energy "dQ" [ \f$ W \cdot m^{-2} \f$ ]
 * - Cold  content of entire snowpack "cc" [\f$ J \cdot m^{-2} \f$ ]
 * - Bulk snow temperature "T_s" [K]
 * - Surface exchange layer temperature "T_s_0" [K]
 * - Lower layer temperature "T_s_l" [K]
 * - Was an iteration error hit "dead". Diagnostic, don't use. [-]
 * - Net shortwave radiation at the surface. Diagnostic. "iswr_out" [ \f$ W \cdot m^{-2} \f$ ]
 * - Binary is the snow isothermal "isothermal" [0,1]
 * - Outgoing longwave radiation "ilwr_out" [ \f$ W \cdot m^{-2} \f$ ]
 * - Total snowpack runoff "sum_snowpack_runoff"
 * - Total snowpack sublimation "sum_snowpack_subl"
 * - Total precipitation onto the snowpack. Diagnostic. "sum_snowpack_pcp" [mm]
 * - Total snowpack melt "sum_melt" [mm]
 * - Snowdepth "snowdepthavg" [m]
 * - Snowdepth in the vertical direction (cosine corrected) "snowdepthavg_vert" [m]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "drift_density": 300,
 *       "const_T_g": -4.0,
 *       "use_slope_SWE": true,
 *       "param_snow_compaction": 1,
 *       "max_h2o_vol":0.0001,
 *       "kt_wetsand": 0.08,
 *       "max_active_layer": 0.1,
 *       "z0":0.001,
 *       "z_T":2.6,
 *       "z_u":2.96,
 *       "z_g":0.1,
 *    }
 *
 * .. confval:: drift_density
 * 
 * 	:default:  300 :math:`kg \cdot m^3`
 * 
 * 	Density of snow that is added from drift events
 *
 *  .. confval:: const_T_g
 * 
 * 	:default:  -4.0 :math:`{}^\circ C`
 * 
 * 	Constant ground temperature
 *
 *  .. confval:: use_slope_SWE
 * 
 * 	:default:  true
 * 
 * 	Use slope corrected snowdepth for compaction. I.e., parallel to gravity force.
 *
 *  .. confval:: param_snow_compaction
 * 
 * 	:default:  1
 * 
 * 	Set to 1 to use new Hedrick, et al parameterization for snow compaction. 0 for origianl Snobal.
 *
 *  .. confval:: max_h2o_vol
 * 
 * 	:default: 0.0001
 * 
 * 	Maximum volumetric water content in the snow. Tends to be required to be set quite low.
 *
 *  .. confval:: kt_wetsand
 * 
 * 	:default:  :math:`0.08 W \cdot m^{-2}`
 * 
 * 	Thermal conductivity of wet sand for G flux
 *
 *  .. confval:: max_active_layer
 * 
 * 	:default:  0.1 m
 * 
 * 	Thickness of active layer
 * 
 *  .. confval:: z0
 * 
 * 	:default: 0.001 m
 * 
 * 	Rouighness length
 * 
 *  .. confval:: z_T
 * 
 * 	:default: 2.6 m
 * 
 * 	Height of air temperature
 * 
 *  .. confval:: z_u
 * 
 * 	:default: 2 m
 * 
 * 	Height of wind measurement
 * 
 *  .. confval:: z_g
 * 
 *     :default: 0.1
 *
 *     Depth of ground temperature measurement
 *
 * \endrst
 *
 * **References:**
 * - Marks, D., Domingo, J., Susong, D., Link, T., Garen, D. (1999). A spatially distributed energy balance snowmelt model
 * for application in mountain basins Hydrological Processes  13(12-13), 1935-1959.
 * - Hedrick, A., Marks, D., Havens, S., Robertson, M., Johnson, M., Sandusky, M., Marshall, H., Kormos, P., Bormann, K., Painter, T. (2018).
 * Direct Insertion of NASA Airborne Snow Observatory‐Derived Snow Depth Time Series Into the iSnobal Energy Balance Snow Model Water Resources Research
 * 54(10), 8045-8063. https://dx.doi.org/10.1029/2018wr023190
 * @}
 */
class snobal : public module_base
{
REGISTER_MODULE_HPP(snobal);
public:
    snobal(config_file cfg);

    ~snobal();

    double drift_density; // if we have blowing snow, this is the density of those particles
    double const_T_g; // constant ground temp, degC

    bool use_slope_SWE; // use a slope corrected SWE for compaction eqn

    virtual void run(mesh_elem &face);
    virtual void init(mesh& domain);
    void checkpoint(mesh& domain, netcdf& chkpt);
    void load_checkpoint(mesh& domain, netcdf& chkpt);

};
