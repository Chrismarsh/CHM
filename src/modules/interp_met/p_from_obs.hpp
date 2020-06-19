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

#include "module_base.hpp"
#include <gsl/gsl_fit.h>
#include <vector>
#include <gsl/gsl_combination.h>

/**
 * \ingroup modules met precip
 * @{
 * \class p_from_obs
 * Spatially distributes liquid water precipitation by calculating lapse rates based on observed data on a per-timestep basis
 *
 *
 * **Depends from met:**
 * - Precipitation  "p" [\f$mm \cdot dt^{-1}\f$]
 *
 * Provides:
 * - Lapsed precipitation "p" [\f$mm \cdot dt^{-1}\f$]
 * - Precipitation corrected for triangle slope. If ``"apply_cosine_correction": false``, then no change. "p_no_slope" [\f$mm \cdot dt^{-1}\f$]
 *
 * **Configuration:**
 * - None
 *
 * @}
 */
class p_from_obs : public module_base
{
REGISTER_MODULE_HPP(p_from_obs);
public:
  p_from_obs(config_file cfg);
    ~p_from_obs();
    void run(mesh_elem& face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};
