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

//
// Created by Chris Marsh on 2015-10-13.
//

#pragma once



#include <string>
#include <cstdio>


//http://stackoverflow.com/a/26197300/410074
template <typename... Ts>
std::string str_format (const std::string &fmt, Ts... vs)
{
    unsigned required = std::snprintf(nullptr, 0, fmt.c_str(), vs...) + 1;
        // See comments: the +1 is necessary, while the first parameter
        //               can also be set to nullptr

    char bytes[required];
    std::snprintf(bytes, required, fmt.c_str(), vs...);

    return std::string(bytes);
}