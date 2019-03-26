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
#include <constants/Atmosphere.h>
#include "module_base.hpp"
#include <boost/shared_ptr.hpp>
#include "logger.hpp"
#include <string>
#include <math.h>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>
/**
 * \addtogroup modules
 * @{
 * \class scale_wind_vert
 * \brief Scales wind speed from reference height to defined heights
 *
 * Depends:
 * - U_R [m/s]
 * - Z_R [m] - Height of wind speed measurement/model layer
 *
 */

class scale_wind_vert : public module_base {
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
