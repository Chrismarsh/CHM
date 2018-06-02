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
* \class MS_wind
* \brief Calculates wind speed and direction following Mason Sykes
*
* Calculates windspeeds using the Mason Sykes wind speed from Essery 1999, using lookup map from DBSM.
 * Calculates the speedup in the mean wind direction defined by all the stations
* Depends:
* - Wind at reference height "U_R" [m/s]
* - Direction at reference height 'vw_dir' [degrees]
*
* Provides:
* - Wind "U_R" [m/s] at reference height
* - Wind direction 'vw_dir' [degrees]
*/
class MS_wind : public module_base
{
public:
    MS_wind(config_file cfg);
    ~MS_wind();
    virtual void run(mesh domain);
    virtual void init(mesh domain);
    double ys;
    double yc;
    class data : public face_info
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
    bool use_ryan_dir;
};

/**
@}
*/