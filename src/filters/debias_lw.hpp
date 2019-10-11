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
#include <constants/Atmosphere.h>


/**
 * \addtogroup filters
 * @{
 * \class debias_lw
 * \brief Debias GEM incoming lw 
 *
 *
 * Depends:
 *
 *
 * Provides:
 *
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
    void init(std::shared_ptr<station>& station);
    void process(std::shared_ptr<station>& station);
};
