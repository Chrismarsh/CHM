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
 * \ingroup modules met wind
 * @{
 * \class uniform_wind
 * Spatially interpolates wind speed and direction without any terrain speed-up modification
 *
 * **Depends from met:**
 * - Wind at reference height "U_R" [ \f$ m \cdot s^{-1}\f$ ]
 * - Direction at reference height "vw_dir" [degrees]
 *
 * **Provides:**
 * - Wind speed at reference height "U_R" [ \f$ m \cdot s^{-1}\f$ ]
 * - Wind direction 'vw_dir' at reference height [degrees]
 *
 */
class uniform_wind : public module_base
{
REGISTER_MODULE_HPP(uniform_wind);
public:
    uniform_wind(config_file cfg);
    ~uniform_wind();
    virtual void run(mesh& domain);
    virtual void init(mesh& domain);
    class lwinddata : public face_info
    {
    public:
        double curvature;
        interpolation interp;
        double corrected_theta;
        double W;
    };
};

/**
@}
*/
