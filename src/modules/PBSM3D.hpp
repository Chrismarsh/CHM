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
//#define BOOST_MATH_INSTRUMENT
#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "interpolation.hpp"

#include "math/coordinates.hpp"

#include <constants/PhysConst.h>
#include "constants/Atmosphere.h"

#include <meteoio/MeteoIO.h>

#include <cmath>
#include <vector>
#include <gsl/gsl_sf_lambert.h>


#ifdef _OPENMP
#include <omp.h>
#else
// we need to define these and have them return constant values under a no-omp situation
inline int omp_get_thread_num() { return 0;}
inline int omp_get_max_threads() { return 1;}

#undef VIENNACL_WITH_OPENMP

#endif


#include <viennacl/linalg/gmres.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/ilu.hpp>
#include <viennacl/linalg/cg.hpp>



#include <armadillo>

#include <cstdlib>
#include <string>
//#define _USE_MATH_DEFINES
//#include <math.h>

#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/roots.hpp>

/**
* \addtogroup modules
* @{
* PBSM3D blowing snow implements a 3D scalar transport model using some PBSM formulations
 * scalar transport using snow concentrations is outlined here
 * Lehning, M., H. Löwe, M. Ryser, and N. Raderschall (2008), Inhomogeneous precipitation distribution and snow transport in steep terrain, Water Resour. Res., 44(7), 1–19, doi:10.1029/2007WR006545.
 * saltation concentration is given in
 * Pomeroy, J. W., and D. M. Gray (1990), Saltation of snow, Water Resour. Res., 26(7), 1583–1594.
 * but clearly given in
 * Pomeroy, J. W., and D. Male (1992), Steady-state suspension of snow, J. Hydrol., 136(1–4), 275–301, doi:10.1016/0022-1694(92)90015-N. [online] Available from: http://linkinghub.elsevier.com/retrieve/pii/002216949290015N
 * deposition flux shown as divergence taken from
 * Liston, G., and M. Sturm (1998), A snow-transport model for complex terrain, J. Glaciol., 44(148), 498–516.
 * sublimation update taken from
 * Pomeroy, J. W., and L. Li (2000), Prairie and arctic areal snow cover mass balance using a blowing snow model, J. Geophys. Res., 105(D21), 26619–26634, doi:10.1029/2000JD900149. [online] Available from: http://www.agu.org/pubs/crossref/2000/2000JD900149.shtml
 * Depends:
* - Air temp (t), [C]
 * - Relative humidity (RH), [%]
 * - SWE [kg/m^2], can come before a snowmodel though.
*
* Provides:
* - mass_drift (kg/m^2) Positive for mass deposition, negative for mass removal.
*/
class PBSM3D : public module_base
{
REGISTER_MODULE_HPP(PBSM3D);
public:
    PBSM3D(config_file cfg);
    ~PBSM3D();
    void run(mesh domain);
    void init(mesh domain);

    double nLayer;
    double susp_depth;
    double v_edge_height;

    // Beta * K, this is beta and scales the eddy diffusivity
    double snow_diffusion_const ;
    double l__max; // vertical mixing length (m)
    bool rouault_diffusion_coeff; //use the spatially variable diffusivity coefficient of Rouault 1991

    bool do_fixed_settling; // should we have a constant settling velocity? 
                            // true: constant settling velocity = settling_velocity (see below)
                            // false: use the parameterization of Pomeroy et al. (1993) and Pomeroy and Gray (1995):
                            //        In this case, the settling velocity decreases with height above the snow surface 
                            //        due to a decrease with height in the mean particle size. 
    double settling_velocity; // Variable used if do_fixed_settling = true
    double n_non_edge_tri;
    double eps; //lapacian smoothing epilson.
    bool limit_mass; // do not saltate more snow than what exists in a cell.
    bool do_sublimation; // should we have a sink sublimation term?
    bool do_lateral_diff; // should have lateral diffusion
    bool enable_veg; // should we consider vegetation ?

    bool use_PomLi_probability; // Use areal Pomeroy Li 2000 probability function.
    bool use_exp_fetch; // Enable the exp Liston 2006 fetch
    bool use_tanh_fetch; // Enable the tanh Pomeroy and Male 1986 fetch


    bool iterative_subl; // if True, enables the iterative sublimation calculation as per Pomeroy and Li 2000
    bool use_R94_lambda; //use the ﻿Raupach 1990 lambda expression using LAI/2 instead of pomeroy stalk density

    double N; //vegetation number density
    double dv; //stalk diameter

    // this is the suspension transport matrix
    double nnz; //number none zero
    viennacl::compressed_matrix<vcl_scalar_type>  vl_C;
    viennacl::vector<vcl_scalar_type> b;

    //this is the drift matrix
    double nnz_drift; //number none zero
    viennacl::compressed_matrix<vcl_scalar_type>  vl_A;
    viennacl::vector<vcl_scalar_type> bb;

    double debug_output;
    double cutoff; // cutoff veg-snow diff (m) that we inhibit saltation entirely
    // don't allow transport if below this threshold.
    // This gives models like snobal a chance to build up their snowpack and avoid convergence issues with thin snowcovers
    double min_mass_for_trans;
    class data : public face_info
    {
    public:
        //edge unit normals
        arma::vec m[5];

        //prism areas
        double A[5];

        //face neighbours
        bool face_neigh[3];

        std::vector<double> u_z_susp; //suspension layer windspeeds
        size_t cell_id;

        double CanopyHeight;
        double LAI;

        //saltation height
        double hs;

        bool is_edge;

        //used to flag the large vegetation areas or other via landcover types to not do any saltation at this point.
        bool saltation;

        double z0;

        double sum_drift;

        double sum_subl;
    };


};

/**
@}
*/
