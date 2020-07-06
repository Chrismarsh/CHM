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
 * \ingroup modules ag
 * @{
 * \class crop_rotation
 *
 * Basic crop rotation. Changes between two different landcovers on even/odd years.
 * On even years, it sets the "crop" parameter to the parameter "annual_crop_inventory_1", and on odd years
 * it sets it to "annual_crop_inventory_2".
 * Mostly a proof-of-concept.
 *
 * **Depends:**
 * - None
 *
 * **Provides:**
 * - None
 *
 * **Configuration:**
 * - None
 *
 * **Parameters:**
 * - Year even crop "annual_crop_inventory_1"
 * - Year odd crop "annual_crop_inventory_2"
 * @}
 */
class crop_rotation : public module_base
{
REGISTER_MODULE_HPP(crop_rotation)
public:
    crop_rotation(config_file cfg);

    ~crop_rotation();

    void run(mesh_elem &face);

};
