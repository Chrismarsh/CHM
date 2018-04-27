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

/**
 * Impliments the fetch algorithm of
 * Lapen, D. R., and L. W. Martz (1993), The measurement of two simple topographic indices of wind sheltering-exposure from raster digital elevation models, Comput. Geosci., 19(6), 769–779, doi:10.1016/0098-3004(93)90049-B.
 */
class fetchr : public module_base
{
public:
    fetchr(config_file cfg);

    ~fetchr();

    virtual void run(mesh_elem& face);

//number of steps along the search vector to check for a higher point
    int steps;
    //max distance to search
    double max_distance;
    double h_IBL; // IBL depth to reestablish steady state (5m to fit blowins snow assumption)
    //size of the step to take
    double size_of_step;

    bool incl_veg;

    //Obstacle heigh increment (m/m)
    //0.06 m/m corresponds to prarie shelter belts
    double I;

};



