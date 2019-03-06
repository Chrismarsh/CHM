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
#include "math/coordinates.hpp"
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <cmath>


#include <viennacl/linalg/gmres.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>

/**
* \addtogroup modules
* @{
* \class Winstral parameter
* \brief Compute the Winstral parameter that can be used to redistribute
*     wind and precipitation in complex terrain to improve simulation of snow 
*     cover variability. 
*
* Depends:
* - Direction at reference height 'vw_dir' [degrees]
* - Snow depth [m] (optional)
*
* Provides:
* - Sx, Winstral paramter [degrees]
*/
class Winstral_parameters : public module_base
{
REGISTER_MODULE_HPP(Winstral_parameters);
public:
    Winstral_parameters(config_file cfg);
    ~Winstral_parameters();


    virtual void run(mesh domain);

    //number of steps along the search vector to check for a higher point
    int steps;
    //max distance to search [m]
    double dmax;
    //size of steps along the search vector  [m]
    double size_of_step;

    // Option to compute the elevation of the point considered to compute Sx
    //     False: the elevation of the center of the triangle is used
    //     True (Recommended): the actual elevation of the point in the plan of the triangle is used
    //        - Improve estimation of Sx around ridges
    //        - Provide smoother estimation of Sx
    bool use_subgridz;

    // height parameter [m] to account for instrument height
    // and/or to reduce the impact of small proximital terrain perturbations
    double height_param;

    // Angular width of the upwind window centered on the wind direction
    // used to compute the average Sx:
    //       - use 0 deg for a Sx just computed in the upwind direction
    //       - use 30 deg as the default value in Winstral et al (2002)
    double angular_window;
    // Angular Increment in the upwind window
    double delta_angle;
    // number of increments in the upwind window
    double nangle;

    // Include local vegetation height when computing Sx
    // Default: False
    bool incl_veg;
    // Include snow deph when computing Sx
    // Improve estimation of Sx when snow is accumulating during the snow season
    bool incl_snw;

    // Calculates the Sx parameter
    double Sx(const mesh &domain, mesh_elem face) const;
};
