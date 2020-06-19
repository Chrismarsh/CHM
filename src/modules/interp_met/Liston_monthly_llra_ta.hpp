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
 * \ingroup modules tair met
 * @{
 * \class Liston_monthly_llra_ta
 *
 * Constant monthly linear lapse rate adjustment for air temperature. Uses Liston and Elder (2006).
 *
 * **Depends from met:**
 * - Air temperature "t" [\f$ {}^\circ C \f$ ]
 *
 * **Provides:**
 * - Air temperature "t" [\f$ {}^\circ C \f$ ]
 * - Lapse rate "t_lapse_rate" [\f$ {}^\circ C \cdot m^{-1} \f$ ]
 *
 * Reference:
 * Liston, Glen E., and Kelly Elder. 2006. A meteorological distribution system for high-resolution terrestrial modeling (MicroMet).
 * Journal of hydrometeorology 7: 217-234
 * @}
 */
class Liston_monthly_llra_ta : public module_base
{
REGISTER_MODULE_HPP(Liston_monthly_llra_ta);
public:
    Liston_monthly_llra_ta(config_file cfg);
    ~Liston_monthly_llra_ta();
    virtual void run(mesh_elem& face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};
