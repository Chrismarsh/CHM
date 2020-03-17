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
    extern void allocate();
    void fsm2_timestep(
    // Driving variables
        float* dt, float* elev, float* zT, float* zU, float* LW, float* Ps, float* Qa, float* Rf, float* Sdif, float* Sdir, float* Sf, float* Ta, float* Ua,

    // Vegetation characteristics
        float* alb0, float* hveg, float* VAI,

    // State variables
        float* albs, float* Tsrf, float* Dsnw, float* Nsnow, float* Qcan, float* Rgrn, float* Sice,
        float* Sliq, float* Sveg, float* Tcan, float* Tsnow, float* Tsoil, float* Tveg, float* Vsmc,

    // Diagnostics
        float* H, float* LE, float* LWout, float* LWsub, float* Melt, float* Roff, float* snd, float* Svg, float* SWE, float* SWout,
        float* SWsub, float* Usub);

    /// module LAYERS
    // These are currently spatially constant and CANNOT be changed on a per triangle basis
    extern  int __layers_MOD_ncnpy;
    extern int __layers_MOD_nsmax;
    extern int __layers_MOD_nsoil;

    extern float __parameters_MOD_rgr0;

    extern float __constants_MOD_e0;
    extern float __constants_MOD_eps;

}

/**
 * \addtogroup modules
 * @{
 * \class FSM
 * \brief FSM Essery, 2019
 *
 *
 * Depends:
 *
 *
 * Provides:
 *
 */
class FSM : public module_base
{
  REGISTER_MODULE_HPP(FSM);
  private:

    struct data : public face_info
    {
        struct
        {
            // Vegetation characteristics
            float alb0 = 0.2;
            float hveg = 0;
            float VAI = 0;
        } veg;

        struct
        {
            // State variables
            float albs = 0.8;
            float Tsrf = 285;
            float Dsnw[3] = {0, 0, 0};
            float Nsnow = 0;
            float Qcan[2] = {0, 0};
            float Rgrn[3] = {__parameters_MOD_rgr0, __parameters_MOD_rgr0, __parameters_MOD_rgr0};
            float Sice[3] = {0, 0, 0};

            float Sliq[3] = {0, 0, 0};
            float Sveg[2] = {0, 0};
            float Tcan[2] = {285, 285};
            float Tsnow[3] = {273, 273, 273};
            float Tsoil[4] = {285, 285, 285, 285};
            float Tveg[2] = {285, 285};

            float Vsat = 0.27;
            float Vsmc[4] = {(float)0.5 * Vsat, (float)0.5 * Vsat, (float)0.5 * Vsat, (float)0.5 * Vsat};

        } state;

        struct
        {
            // Diagnostics
            float H = -9999;
            float LE = -9999;
            float LWout = -9999;
            float LWsub = -9999;
            float Melt = -9999;
            float Roff = -9999;

            float snd = -9999;
            float Svg = -9999;
            float SWE = -9999;
            float SWout = -9999;

            float SWsub = -9999;
            float Usub = -9999;
        } diag;
    };

  public:
    FSM(config_file cfg);
    ~FSM();
    virtual void run(mesh_elem& face);
    virtual void init(mesh& domain);
};

/**
@}
*/
