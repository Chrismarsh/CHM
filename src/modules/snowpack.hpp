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

#include <memory>

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"



#include <meteoio/MeteoIO.h>
#include <physics/Atmosphere.h>
#include <snowpack/libsnowpack.h>

#include <string>

/**
 * \ingroup modules snow
 * @{
 * \class Lehning_snowpack
 *
 * SNOWPACK (Bartelt and Lehning, 2002) is a multilayer finite-element energy balance snow model originally developed for avalanche hazard forecasting.
 * It has the greatest computational burden of all snowpack models in CHM. This version of SNOWPACK has been modified to allow unlimited removal of snow
 * layers by the blowing snow and avalanche routines.
 *
 * \rst
 * .. note::
 *    Currently SNOWPACK is setup to be run an external albedo model and should be changed to default back to the SNOWPACK one.
 * \endrst
 *
 * **Depends:**
 * - Incoming shortwave radiation, all beams "iswr" [ \f$ W \cdot m^{-2} \f$ ]
 * - Incomging longwave radiation "ilwr" [ \f$ W \cdot m^{-2} \f$ ]
 * - Relative humidity "rh" [%]
 * - Air temperature "t" [\f$ {}^\circ C \f$]
 * - Windspeed at 2m "U_2m_above_srf" [\f$ m \cdot s^{-1}    \f$]
 * - Precipitation "p" [\f$ mm \cdot dt^{-1} \f$]
 * - Precipitation rain fraction "frac_precip_rain" [-]
 * - Snow albedo "snow_albedo" [-]
 *
 * **Optional:**
 * - Optionally, depend on the ``_subcanopy`` variants of the above
 * - Blowing snow erosion/deposition mass "drift_mass" [mm]
 * - Ground temperature "T_g" [\f$  {}^\circ C   \f$]
 *
 * **Provides:**
 * - Change in internal energy "dQ" [ \f$ W \cdot m^{-2} \f$ ]
 * - Snow water equivalent "swe" [mm]
 * - Bulk snow temperature "T_s" [\f$ {}^\circ C \f$]
 * - Surface exchange layer temperature "T_s_0"
 * - Number of discretization nodes "n_nodes"
 * - Number of discretization layers "n_elem"
 * - Snowdepth "snowdepthavg" [m]
 * - Sensible heat flux "H" [ \f$ W \cdot m^{-2} \f$ ]
 * - Latent heat flux "E" [ \f$ W \cdot m^{-2} \f$ ]
 * - Ground heat flux "G" [ \f$ W \cdot m^{-2} \f$ ]
 * - Outgoing longwave radiation "ilwr_out" [ \f$ W \cdot m^{-2} \f$ ]
 * - Reflected shortwave radiation "iswr_out" [ \f$ W \cdot m^{-2} \f$ ]
 * - Net allwave radiation "R_n" [ \f$ W \cdot m^{-2} \f$ ]
 * - Snowpack runoff "runoff" [mm]
 * - Snow mass removed "mass_snowpack_removed" [mm]
 * - Total snowpack runoff "sum_runoff"
 * - Total surface sublimation loss "sum_subl"
 * - Sublimation loss "sublimation"
 * - Evaporation loss "evap"
 * - "MS_SWE"
 * - "MS_WATER"
 * - "MS_TOTALMASS"
 * - "MS_SOIL_RUNOFF"
 * - Top layer Snow sphericity "snow_sphericity" [-]
 * - Top layer Snow grain size "snow_grain_size" [mm]
 * - Top layer Fraction of ice content "frac_ice_content" [-]
 * - Top layer Liquid water content "Sliq" [-]
 * - Snow surface temperature "Tsnow" [\f$ {}^\circ C \f$]
 * - Top layer Snow density "snow_density" [\f$ kg \cdot m^{-3} \f$]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    "Lehning_snowpack":
 *    {
 *      "sno":
 *      {
 *         "SoilAlbedo": 0.09,
 *         "BareSoil_z0": 0.2,
 *         "WindScalingFactor": 1,
 *         "TimeCountDeltaHS": 0.0
 *      },
 *      "GROOMING_WEEK_START": 40,
 *      "GROOMING_WEEK_END": 17,
 *      "GROOMING_HOUR": 21,
 *      "GROOMING_DEPTH_START": 0.4,
 *      "GROOMING_DEPTH_IMPACT": 0.4
 *    }
 * \endrst
 * @}
 */
class Lehning_snowpack : public module_base
{
REGISTER_MODULE_HPP(Lehning_snowpack);
public:
    Lehning_snowpack(config_file cfg);

    ~Lehning_snowpack();

    virtual void run(mesh_elem &face);

    virtual void init(mesh& domain);

    virtual void checkpoint(mesh& domain, netcdf& chkpt);
    virtual void load_checkpoint(mesh& domain, netcdf& chkpt);

    struct data : public face_info
    {
        //main snowpack model
        std::shared_ptr<Snowpack> sp;

        /*
         * This is the PRIMARY data structure of the SNOWPACK program \n
         * It is used extensively not only during the finite element solution but also to control
         */
        std::shared_ptr<SnowStation> Xdata;

        std::shared_ptr<SnowpackConfig> Spackconfig;
        std::shared_ptr<Meteo> meteo;
        std::shared_ptr<Stability> stability;
        mio::Config config;
        double cum_precip;

        double sum_subl;
        bool is_grooming; // Flag indicating if snow grooming is enabled for this face
    };

    double sn_dt; // calculation step length
    double const_T_g; // constant ground temp, degC

};
