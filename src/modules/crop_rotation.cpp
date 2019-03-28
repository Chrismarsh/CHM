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

#include "crop_rotation.hpp"
REGISTER_MODULE_CPP(crop_rotation);

crop_rotation::crop_rotation(config_file cfg)
        : module_base("crop_rotation", parallel::data, cfg)
{
    provides_parameter("crop"_s);
}

crop_rotation::~crop_rotation()
{

}

void crop_rotation::run(mesh_elem& face)
{

    int year = global_param->year();

    if(year == 2010)
        face->parameter("crop"_s) = face->parameter("annual_crop_inventory_2010"_s);
    else if( year == 2011)
        face->parameter("crop"_s) = face->parameter("annual_crop_inventory_2011"_s);

}
