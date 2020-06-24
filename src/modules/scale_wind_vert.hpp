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

#include <boost/shared_ptr.hpp>
#include <physics/Atmosphere.h>
#include "module_base.hpp"
#include <boost/shared_ptr.hpp>
#include "logger.hpp"
#include <string>
#include <math.h>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>
/**
 * \ingroup modules wind
 * @{
 * \class scale_wind_vert
 * Scales wind speed from reference height to defined heights. Scales the wind vertically using a neutral stability
 * log relationship. If snow depth is present, the height is taken into account. A fixed z0 of 0.01m is currently used.
 * If a canopy is present then an exp sub-canopy scaling is used.
 *
 * **Depends:**
 * - Wind speed @reference height "U_R" [ \f$ m \cdot s^{-1} \f$ ]
 *
 * **Optional:**
 * - Snow depth to take into account for the scaling "snowdepthavg" [m]
 *
 * **Provides:**
 * - Wind speed @2m "U_2m_above_srf" [ \f$ m \cdot s^{-1} \f$ ]
 *
 * **Parameters:**
 * - If ``ignore_canopy=False``, then Leaf Area Index "LAI" [-]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "ignore_canopy": false
 *    }
 *
 * .. confval:: ignore_canopy
 *
 *    :default: false
 *
 *    Ignores the vegetation canopy
 *
 * \endrst
 *
 * @}
 */
class scale_wind_vert : public module_base
{
REGISTER_MODULE_HPP(scale_wind_vert);
public:
    scale_wind_vert(config_file cfg);

    ~scale_wind_vert();
    virtual void init(mesh& domain);

    //this module can swap between a domain parallel and a data parallel state
    ///domain parallel allows for blending through vegetation to avoid sharp gradietns that can complicate blowing snow, &c.
    virtual void run(mesh& domain);
    virtual void run(mesh_elem &face);

    //scales the windspeed at a single triangle
    void point_scale(mesh_elem &face);

    bool ignore_canopy;
    //virtual void init(mesh& domain);
    struct d: public face_info
    {
        double temp_u;
        interpolation interp;
    };
};
