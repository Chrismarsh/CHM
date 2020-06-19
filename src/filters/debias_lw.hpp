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

#include "filter_base.hpp"
#include <physics/Atmosphere.h>


/**
 * \ingroup filters lw
 *
 * @{
 * \class debias_lw
 * Debias an incoming longwave radiation using an additive constant factor.
 *
 *
 * **Requires:**
 * - Incoming longwave radiation \f$[W \cdot m^{-2}\f$]
 *
 * **Modifies:**
 * - Input longwave
 *
 * **Configuration keys:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "variable": "lw",
 *       "factor": 3.5
 *    }
 *
 * .. confval:: variable
 *
 *    :type: string
 *
 *    Name of the variable to modify that coincides with the input met file.
 *
 * .. confval:: factor
 *
 *    :type: double
 *    :units: :math:`[W \cdot m^{-2}]`
 *
 *    The amount to add to the input longwave by: `lw = lw + factor`
 *
 * \endrst
 * @}
 */
class debias_lw : public filter_base
{
REGISTER_FILTER_HPP(debias_lw);
private:
    std::string var;
    double fac;
public:
    debias_lw(config_file cfg);
    ~debias_lw();
    void init();
    void process(std::shared_ptr<station>& station);
};

