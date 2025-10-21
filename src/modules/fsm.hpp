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
#include <string>

extern "C"
{
    void fsm_allocate();
    void fsm_layers_config(float fvg1, float zsub, int ncnpy, int nsmax, int nsoil);
    void fsm_soilprops_config(float b, float hcap_soil, float hcon_soil, float sathh, float vcrit, float vsat);
    float fsm_constants_e0();
    float fsm_constants_eps();
    float fsm_parameters_rgr0();
    void fsm2_timestep(
    // Driving variables
        float* dt, float* elev, float* zT, float* zU,
        float* LW, float* Ps, float* Qa, float* Rf, float* Sdif, float* Sdir, float* Sf, float* Ta, float* trans, float* Ua,

    // Vegetation characteristics
        float* alb0, float* hveg, float* VAI,

        float* rhod,

    // State variables
        float* albs, float* Tsrf, float* Dsnw, int* Nsnow, float* Qcan, float* Rgrn, float* Sice,
        float* Sliq, float* Sveg, float* Tcan, float* Tsnow, float* Tsoil, float* Tveg, float* Vsmc,

    // Diagnostics
        float* H, float* LE, float* LWout, float* LWsub, float* Melt, float* Roff, float* snd, float* snw, float* subl, float* svg,
        float* SWout, float* SWsub, float* Usub, float*  Wflx
        );
}

/**
 * \ingroup modules snow
 * @{
 * \class FSM
 * Flexible Snow Model (FSM) 2.0
 *
 * "The Flexible Snow Model (FSM2) is a multi-physics energy balance model of snow accumulation and melt,
 * extending the Factorial Snow Model (Essery, 2015) with additional physics, driving and output options."
 *
 * This version of FSM has been customized to have the sophisticated process parametrizations selected in every case,
 * except atmospheric stability corrections.
 *
 *  * \rst
 * .. warning::
 *    Vegetation characteristics remain a TODO
 *
 * \endrst
 *
 * **Depends:**
 * - Solar elevation "solar_el" [degrees]
 * - Incoming longwave radiation "ilwr" \f$[W \cdot m^{-2}\f$]
 * - Relative Humidy "rh" [%]
 * - Air temperature "t" [ \f$ {}^\circ C \f$]
 * - Precipitation snow "p_snow" [\f$mm \cdot dt^{-1}\f$]
 * - Precipitation rain "p_rain" [\f$mm \cdot dt^{-1}\f$]
 * - Wind speed 2 m above surface "U_2m_above_srf" [ \f$ m \cdot s^{-1}\f$ ]
 * - Incoming shortwave radiation, direct beam "iswr_direct" \f$[W \cdot m^{-2}\f$]
 * - Incoming shortwave radiation, diffuse beam "iswr_diffuse" \f$[W \cdot m^{-2}\f$]
 *
 * **Provides:**
 * - Sensible heat flux "H" \f$[W \cdot m^{-2}]\f$
 * - Latent heat flux "E" \f$[W \cdot m^{-2}]\f$
 * - Outgoing longwave radiation "ilwr_out" (alias "LWout") \f$[W \cdot m^{-2}]\f$
 * - Melt rate "melt_rate" \f$[mm \cdot dt^{-1}]\f$
 * - Surface runoff/drainage, includes rain on snow percolation "roff" \f$[mm \cdot dt^{-1}]\f$
 * - Mean snow depth (slope-normal) "snowdepthavg" [m]
 * - Mean snow depth (vertical projection) "snowdepthavg_vert" [m]
 * - Snow water equivalent "swe" [mm]
 * - Instantaneous sublimation "subl" \f$[mm \cdot dt^{-1}]\f$
 * - Cumulative sublimation "sum_snowpack_subl" [mm]
 * - Snow albedo "snow_albedo" [-]
 * - Number of snow layers "Nsnow" [-]
 * - Snow liquid water content layer 0 "Sliq[0]" [mm]
 * - Snow liquid water content layer 1 "Sliq[1]" [mm]
 * - Snow liquid water content layer 2 "Sliq[2]" [mm]
 * - Snow temperature layer 0 "Tsnow[0]" [K]
 * - Snow temperature layer 1 "Tsnow[1]" [K]
 * - Snow temperature layer 2 "Tsnow[2]" [K]
 * - Soil temperature layer 0 "Tsoil[0]" [K]
 * - Soil temperature layer 1 "Tsoil[1]" [K]
 * - Soil temperature layer 2 "Tsoil[2]" [K]
 * - Soil temperature layer 3 "Tsoil[3]" [K]
 *
 * **Optional:**
 *
 * Sub-canopy forcing:
 * - Subcanopy incoming shortwave radiation "iswr_subcanopy" \f$[W \cdot m^{-2}\f$]
 * - Subcanopy relative humidity "rh_subcanopy" [%]
 * - Subcanopy air temperatue "ta_subcanopy" [ \f$ {}^\circ C \f$]
 * - Subcanopy incoming longwave radidation "ilwr_subcanopy" \f$[W \cdot m^{-2}\f$]
 *
 *
 * **References:**
 * - https://github.com/RichardEssery/FSM2
 * - Essery, R. (2015). A factorial snowpack model (FSM 1.0) Geoscientific Model Development  8(12), 3867 3876.
 * https://dx.doi.org/10.5194/gmd-8-3867-2015
 *
 * @}
 */
class FSM : public module_base
{
  REGISTER_MODULE_HPP(FSM);
  private:
    static constexpr int SNOW_LAYER_COUNT = 6;
    static constexpr int SOIL_LAYER_COUNT = 4;
    static constexpr int CANOPY_LAYER_COUNT = 2;

    struct data : public face_info
    {
        struct
        {
            // Vegetation characteristics
            float alb0 = 0.2;
            float vegh = 0;
            float VAI = 0;
        } veg;

        struct
        {
            // State variables
            float albs = 0.8;
            float Tsrf = 263.0; // cold soils
            float Dsnw[SNOW_LAYER_COUNT] = {0, 0, 0, 0, 0, 0}; //Snow layer thicknesses (m)
            int Nsnow = 0; //Number of snow layers
            float Qcan[CANOPY_LAYER_COUNT] = {0, 0}; // Canopy air space humidities
            float Rgrn[SNOW_LAYER_COUNT] = {0, 0, 0, 0, 0, 0}; //Snow layer grain radii (m)
            float Sice[SNOW_LAYER_COUNT] = {0, 0, 0, 0, 0, 0}; // Ice content of snow layers (kg/m^2)

            float Sliq[SNOW_LAYER_COUNT] = {0, 0, 0, 0, 0, 0}; // Liquid content of snow layers (kg/m^2)
            float Sveg[CANOPY_LAYER_COUNT] = {0, 0}; //Snow mass on vegetation layers (kg/m^2)
            float Tcan[CANOPY_LAYER_COUNT] = {285, 285}; //Canopy air space temperatures (K)
            float Tsnow[SNOW_LAYER_COUNT] = {263, 263, 263, 263, 263, 263}; //Snow layer temperatures (K)
            float Tsoil[SOIL_LAYER_COUNT] = {263, 263.1, 263.2, 263.3}; //Soil layer temperatures (K); Cold soils
            float Tveg[CANOPY_LAYER_COUNT] = {285, 285}; // Vegetation layer temperatures (K)

            float Vsat = 0.27;

            // Volumetric moisture content of soil layers
            float Vsmc[SOIL_LAYER_COUNT] = {(float)0.5 * Vsat, (float)0.5 * Vsat, (float)0.5 * Vsat, (float)0.5 * Vsat};

        } state;

        struct
        {
            // Diagnostics
            float H = -9999; // Sensible heat flux to the atmosphere (W/m^2)
            float LE = -9999; // Latent heat flux to the atmosphere (W/m^2)
            float LWout = -9999; // Outgoing LW radiation (W/m^2)
            float LWsub = -9999; // Subcanopy downward LW radiation (W/m^2)
            float Melt = -9999; // Surface melt rate (kg/m^2/s)
            float Roff = -9999; // Runoff from snow (kg/m^2/s)

            float snd = -9999; // Snow depth (m)
            float snw = -9999; // Total snow mass on ground (kg/m^2)
            float subl = -9999; // Sublimation rate (kg/m^2/s)
            float svg = -9999; // Total snow mass on vegetation (kg/m^2)

            float SWout = -9999; // Outgoing SW radiation (W/m^2)
            float SWsub = -9999; // Subcanopy downward SW radiation (W/m^2)
            float Usub = -9999; // Subcanopy wind speed (m/s)
            float Wflx[SNOW_LAYER_COUNT] = {-9999, -9999, -9999, -9999, -9999, -9999}; // Water flux into snow layer (kg/m^2/s)

            float sum_snowpack_subl = -9999; // cumulative sublimation (kg/m^2)
        } diag;
    };

    float constants_e0_ = 0.0f;
    float constants_eps_ = 0.0f;
    float parameters_rgr0_ = 0.0f;

  public:
    FSM(config_file cfg);
    ~FSM();
    virtual void run(mesh_elem& face);
    virtual void init(mesh& domain);
    void checkpoint(mesh& domain, netcdf& chkpt);
    void load_checkpoint(mesh& domain, netcdf& chkpt);
};
