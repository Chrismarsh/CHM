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
#include <math.h>


/**
 * \addtogroup filters
 * @{
 * \class goodison_undercatch
 * \brief Computes undercatch correction
 *
 * Undercatch correction for a Nipher shielded guage via Goodison 1998 for solid precipitation
 *
 * Depends:
 * - p [mm]
 * - u [m/s]
 *
 * References:
 * - ﻿Goodison, B. E. (1998), WMO Solid Solid Precipitiation Measurement Intercomparison. https://globalcryospherewatch.org/bestpractices/docs/WMOtd872.pdf
 * @}
 * **/
class goodison_undercatch : public filter_base
{
REGISTER_FILTER_HPP(goodison_undercatch);
private:
    std::string var;
public:
    goodison_undercatch(config_file cfg);
    ~goodison_undercatch();
    void init();
    void process(std::shared_ptr<station>& station);
};
