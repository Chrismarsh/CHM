//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

/**
 * \ingroup modules met rh
 * @{
 * \class rh_no_lapse
 * Spatially interpolates relative humidity without any lapse adjustments. Bounds RH on [10%, 100%].
 *
 * **Depends from met:**
 * - Relative Humidity "rh" [%]
 *
 * **Provides:**
 * - Relative Humidity "rh" [%]
 *
 * **Configuration keys:**
 * - None
 *
 * @}
 */
class rh_no_lapse : public module_base
{
REGISTER_MODULE_HPP(rh_no_lapse);
public:
    rh_no_lapse(config_file cfg);

    ~rh_no_lapse();

    virtual void run(mesh_elem &face);
    virtual void init(mesh& domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};



