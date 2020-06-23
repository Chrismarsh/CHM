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
#include "TPSpline.hpp"
#include <meteoio/MeteoIO.h>

/**
 * \ingroup modules iswr
 * @{
 * \class Iqbal_iswr
 *
 * Estimates incoming shortwave direct and diffuse beams without slope correction.
 *
 * **Depends:**
 * - Air temperature "t" [ \f${}^\circ C \f$]
 * - Relative humidity "rh" [%]
 * - Cloud fraction "cloud_frac" [-]
 * - Solar elevation "solar_el" [degrees]
 *
 * **Provides:**
 * - Incoming solar shortwave radiation, direct beam, no slope adjustment "iswr_direct_no_slope" [ \f$ W \cdot m^{-2}\f$]
 * - Incoming solar shortwave radiation, diffuse beam, no slope adjustment "iswr_diffuse_no_slope" [ \f$ W \cdot m^{-2}\f$]
 * - Atmospheric transmittance "atm_trans" [-]
 *
 * **References:**
 * - Code ported from SUMMA
 * - Iqbal, M. (1983). An Introduction to Solar Radiation https://dx.doi.org/10.1016/b978-0-12-373750-2.50015-x
 * @}
 */
class Iqbal_iswr : public module_base
{
REGISTER_MODULE_HPP(Iqbal_iswr);
public:
    Iqbal_iswr(config_file cfg);
    ~Iqbal_iswr();
    void run(mesh_elem& face);
};
