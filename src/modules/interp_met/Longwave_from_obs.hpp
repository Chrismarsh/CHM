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

#include "../logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <cstdlib>
#include <string>

#include <cmath>

#include <math.h>

/**
 * \ingroup modules lw met
 * @{
 * \class Longwave_from_obs
 * Annually constant longwave lapse rate adjustment to longwave from observations.
 * Lapse rate of 2.8 W/m^2 / 100 meters
 *
 * **Depends from met:**
 * - Incoming longwave radiation "Qli"  \f$[W \cdot m^{-2}\f$]
 *
 * **Provides:**
 * - Incoming longwave "ilwr"  \f$[W \cdot m^{-2}\f$]
 *
 * **Configuration keys:**
 * - None
 *
 * **Reference:**
 * Marty, C., Philipona, R., Fröhlich, C., Ohmura, A. (2002).
 * Altitude dependence of surface radiation fluxes and cloud forcing in the alps: results from the alpine surface
 * radiation budget network Theoretical and Applied Climatology  72(3), 137-155. https://dx.doi.org/10.1007/s007040200019
 *
 * @}
*/
class Longwave_from_obs : public module_base
{
REGISTER_MODULE_HPP(Longwave_from_obs);
public:
    Longwave_from_obs(config_file cfg);
    ~Longwave_from_obs();
    virtual void run(mesh_elem& face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};
