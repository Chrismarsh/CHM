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

/**
 * \ingroup modules pointmode met
 * @{
 * \class point_mode
 *
 * Use this module to enable using CHM in point mode. This module does not to any spatial interpolation. Instead, it passes
 * input met through to the per-triangle variables storages. This is suitable for using CHM like a point scale model, such as at a
 * research site. Specifically this sets ``depends_from_met(...)`` inputs such that ``depends(...)`` in other modules can work.
 *
 * \rst
 * .. note::
 *    The configuration allows for fine-tuning what is passed through. So although ``Depends from met`` below lists all of the possible depends,
 *    whatever is set in the configuration will be what is required at runtime.
 *
 * \endrst
 *
 * **Depends from met::**
 * - Air temp "t" [\f$ {}^\circ C \f$]
 * - Relative humidity "rh" [%]
 * - Wind speed "u" [ \f$ m \cdot s^{-1} \f$ ]
 * - Precipitation "p" [mm]
 * - Incoming longwave radiation "Qli" [ \f$ W \cdot m^{-2}\f$]
 * - Incoming shortwave radiation "Qsi" [ \f$ W \cdot m^{-2}\f$]
 * - Incomging shortwave radiation, diffuse beam "iswr_diffuse" [ \f$ W \cdot m^{-2}\f$]
 * - Incomging shortwave radiation, direct beam "iswr_direct" [ \f$ W \cdot m^{-2}\f$]
 * - Wind direction @2m "vw_dir" [degrees]
 * - Ground temperature "T_g" [\f$ {}^\circ C \f$]
 *
 * **Provides:**
 * - Air temp "t" [\f$ {}^\circ C \f$]
 * - Relative humidity "rh" [%]
 * - Wind speed @2m "U_2m_above_srf" [ \f$ m \cdot s^{-1} \f$ ]
 * - Wind speed @reference height "U_R" [ \f$ m \cdot s^{-1} \f$ ]
 * - Precipitation "p" [mm]
 * - Incoming longwave radiation "ilwr" [ \f$ W \cdot m^{-2}\f$]
 * - Incoming shortwave radiation "iswr" [ \f$ W \cdot m^{-2}\f$]
 * - Incomging shortwave radiation, diffuse beam "iswr_diffuse" [ \f$ W \cdot m^{-2}\f$]
 * - Incomging shortwave radiation, direct beam "iswr_direct" [ \f$ W \cdot m^{-2}\f$]
 * - Wind direction @2m "vw_dir" [degrees]
 * - Ground temperature "T_g" [\f$ {}^\circ C \f$]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "t": true,
 *       "rh": true,
 *       "U_R": true,
 *       "U_2m_above_srf": true,
 *       "p": true,
 *       "ilwr": true,
 *       "iswr": true,
 *       "vw_dir": true,
 *       "iswr_diffuse": false,
 *       "iswr_direct": false,
 *       "T_g": false,
 *
 *       "override":
 *       {
 *          "svf":0.8
 *       }
 *
 *   }
 *
 * .. confval:: X
 *    :type: boolean
 *
 *    When ``true``, pass _X_ through from the forcing file.
 *
 * .. confval:: override.svf
 *
 *    Allows for setting the sky view factor at the point to be different than what was calculated from the mesh.
 * \endrst
 * @}
 */
class point_mode : public module_base
{
REGISTER_MODULE_HPP(point_mode);
public:
    point_mode(config_file cfg);

    ~point_mode();

    virtual void run(mesh_elem &face);


    bool t ;
    bool rh ;
    bool U_R ;
    bool p ;
    bool ilwr ;
    bool iswr ;
    bool vw_dir ;
    bool iswr_diffuse ;
    bool iswr_direct ;
    bool U_2m_above_srf ;
    bool T_g;


};
