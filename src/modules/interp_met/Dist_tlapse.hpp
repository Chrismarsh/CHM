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
 * \class Dist_tlapse
 *  Spatially interpolates provided lapse rates from virtual stations
 *
 * **Depends from met:**
 * - Air temperature - "t" [\f$ \circ C ] \f$
 * - Lapse rate - "t_lapse_rate" [\f$ \circ C \cdot m^{-1}] \f$
 *
 * **Provides:**
 * - Air temperature "t" [\f$ \circ C ] \f$
 * - Lapse rate "t_lapse_rate" [\f$ \circ C \cdot m^{-1}] \f$
 *
 * **Configuration keys:**
 * - None
 * @}
 */
class Dist_tlapse : public module_base
{
REGISTER_MODULE_HPP(Dist_tlapse);
public:
    Dist_tlapse(config_file cfg);
    ~Dist_tlapse();
    virtual void run(mesh_elem& face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};

