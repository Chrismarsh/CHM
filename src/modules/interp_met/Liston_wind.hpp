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
* \class Liston_wind
* \brief Calculates wind speed and direction following Liston 2006
*
* Calculates windspeeds using a terrain curvature following Liston 2006.
* Defaults to 300 m distance. Configurable using "distance".
 * Liston, G. E., & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). Journal of Hydrometeorology, 7(2), 217–234. http://doi.org/10.1175/JHM486.1
* Depends:
* - Wind at reference height "U_R" [m/s]
* - Direction at reference height 'vw_dir' [degrees]
*
* Provides:
* - Wind "U_R" [m/s] at reference height
* - Wind direction 'vw_dir' [degrees]
*/
class Liston_wind : public module_base
{
REGISTER_MODULE_HPP(Liston_wind);
public:
    Liston_wind(config_file cfg);
    ~Liston_wind();
    virtual void run(mesh domain);
    virtual void init(mesh domain);
    double ys;
    double yc;
    class lwinddata : public face_info
    {
    public:
        double curvature;
        interpolation interp;
        double corrected_theta;
        double W;
        double temp_u;
        interpolation interp_smoothing;
    };
    double distance;
    double Ww_coeff;
};

/**
@}
*/
