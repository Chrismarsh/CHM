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
#include <meteoio/MeteoIO.h>

/**
 * Dodson, R. and Marks, D.: Daily air temperature interpolated at high spatial resolution over a large mountainous region, Clim. Res., 8(Myers 1994), 1–20, doi:10.3354/cr008001, 1997.
 */
class Dodson_NSA_ta : public module_base
{
public:
    Dodson_NSA_ta(config_file cfg);
    ~Dodson_NSA_ta();
    void run(mesh_elem &face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};


