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
#include "TPSpline.hpp"

// Calculates horizon shadows using an adaptation of Dozier and Frew 1990 to an unstructured mesh
//Dozier, J., & Frew, J. (1990).
// Rapid calculation of terrain parameters for radiation modeling from digital elevation data.
// IEEE Transactions on Geoscience and Remote, 28(5), 963–969.
class fast_shadow : public module_base
{
REGISTER_MODULE_HPP(fast_shadow);
public:
    fast_shadow(config_file cfg);

    ~fast_shadow();

    virtual void run(mesh_elem& face);

//number of steps along the search vector to check for a higher point
    int steps;
    //max distance to search
    double max_distance;

    //size of the step to take
    double size_of_step;

};
