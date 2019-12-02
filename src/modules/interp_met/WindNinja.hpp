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
#include <physics/Atmosphere.h>
#include <cstdlib>
#include <string>

#include <Winstral_parameters.hpp>

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
* \class WindNinja
* \brief Calculates wind speed and direction following the downscaling stategy of Barcons et al. (Wind Energy, 2018)
*
* Calculates windspeeds and direction from GEM input and a library of high-resolution wind field generated with the WindNinja wind flow model
* Depends:
* - Wind at reference height "U_R" [m/s]
* - Direction at reference height 'vw_dir' [degrees]
*
* Provides:
* - Wind "U_R" [m/s] at reference height
* - Wind direction 'vw_dir' [degrees]
*/
class WindNinja : public module_base
{
REGISTER_MODULE_HPP(WindNinja);
public:
    WindNinja(config_file cfg);
    ~WindNinja();
    virtual void run(mesh& domain);
    virtual void init(mesh& domain);
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
        double W_transf;
    };
    double distance;
    int N_windfield; //  Number of wind fields in the library
    bool ninja_average; // Boolean to activate linear interpolation betweem the closest 2 wind fields from the library 
    double H_forc; // Reference height for GEM forcing and WindNinja wind field library
    double Max_spdup;  // Maximal value of crest speedup
    double Min_spdup;  // Minimal value of crest speedup
    int L_avg;   // Size of the square of averaging when computing the speed-up map 
                 // Not used by default and set to -1 if not specified in the config file. 
    bool ninja_recirc; // Boolean to activate wind speed reduction on the leeside of mountainous terrain

    bool compute_Sx; // uses the Sx module to influence the windspeeds so Sx needs to be computed during the windspeed evaluation, instead of a seperate module
    double Sx_crit;    // Critical values of the Winstral parameter to determine the occurence of flow separation.  
    boost::shared_ptr<Winstral_parameters> Sx;
};

/**
@}
*/
