/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <meteoio/MeteoIO.h>
/**
 * \addtogroup modules
 * @{
 * \class Burridge_iswr
 * \brief Computes incoming shortwave radiation using a cloud fraction.
 *
 * Computes incoming direct and diff shortwave radiation using a cloud fraction based on RH at 700mb.
 *
 * Depends:
 * - cloud_frac (-)
 *
 * Provides:
 * - Shortwave all beam "iswr" [W/m^2]
 * - Shortwave direct "iswr_direct" [W/m^2]
 * - Shortwave diffuse "iswr_diffuse" [W/m^2]
 * - Atmospheric transmittance, [0,1] "atm_trans" [-]
 * References:
 * - Burridge, D. M., and A. J. Gadd, 1974: The Meteorological Office operational 10 level numerical weather prediction model (December 1974). U.K. Met. Office Tech. Notes 12 and 48, 57 pp.
 * - Described in Liston, G. E., & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). Journal of Hydrometeorology, 7(2), 217–234. http://doi.org/10.1175/JHM486.1
 */
class Burridge_iswr : public module_base
{
public:
    Burridge_iswr(config_file cfg);
    ~Burridge_iswr();
    void run(mesh_elem& face);

};



