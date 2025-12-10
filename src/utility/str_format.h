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



#include <stdexcept>
#include <string>
#include <cstdio>
#include <vector>

//https://stackoverflow.com/a/26221725
template <typename... Ts>
std::string str_format (const std::string &fmt, Ts... vs)
{
    int size_s = std::snprintf(nullptr, 0, fmt.c_str(), vs ...) +1;
    if (size_s <= 0)
        throw std::runtime_error("Error during string formatting in str_format.h.");

    auto size = static_cast<size_t>(size_s);
    std::vector<char>  buf(size);

    auto result = std::snprintf(buf.data(), size, fmt.c_str(), vs ...);
    if (result < 0 || static_cast<size_t>(result) >= size)
        throw std::runtime_error("Error during string formatting in str_format.h.");

    return std::string(buf.data(), buf.data() + result);
}
