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
 * \addtogroup filters
 * @{
 * \class filter_template_name
 * \brief Template header file for new filter implementations.
 *
 *
 * Depends:
 *
 *
 * Provides:
 * @}
 */
class filter_template : public filter_base
{
REGISTER_FILTER_HPP(filter_template);
private:
    double filter_variable;
public:
    filter_template(config_file cfg);
    ~filter_template();
    void init();
    void process(std::shared_ptr<station>& station);
};
