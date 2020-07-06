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
#include "physics/Atmosphere.h"
#include <cstdlib>
#include <string>
#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

/**
 * \ingroup modules exp evap
 * @{
 * \class PenmanMonteith_evaporation

 * Calculates evapo-transpiration via Penman-Monteith.
 *
 * Not currently maintained.
 *
 * Depends:
 * - Incoming shortwave radaition "iswr" [\f$ W \cdot m^{-2} \f$ ]
 * - Incoming longwave radiation "ilwr" [\f$ W \cdot m^{-2} \f$ ]
 * - Relative humidity "rh" [%]
 * - Windspeed at 2m "U_2m_above_srf" [\f$ m \cdot s^{-1} \f$ ]
 *
 * Provides:
 * - Evapotranspiration "ET" [\f$ mm \cdot dt^{-1} \f$ ]
 *
 * @}
 */
class PenmanMonteith_evaporation : public module_base
{
REGISTER_MODULE_HPP(PenmanMonteith_evaporation);
public:
    PenmanMonteith_evaporation(config_file cfg);
    ~PenmanMonteith_evaporation();
    virtual void run(mesh_elem& face);

};
