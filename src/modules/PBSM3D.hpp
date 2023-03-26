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
#include "interpolation.hpp"
#include "logger.hpp"
#include "module_base.hpp"
#include "triangulation.hpp"

#include "math/coordinates.hpp"
#include "math/LinearAlgebra.hpp"

#include <physics/PhysConst.h>
#include "physics/Atmosphere.h"

#include <meteoio/MeteoIO.h>

#include <cmath>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_lambert.h>
#include <vector>

#include <armadillo>

#include <cstdlib>
#include <string>
//#define _USE_MATH_DEFINES
//#include <math.h>

#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/tuple.hpp>


/**
 * \ingroup modules snow
 * @{
 * \class PBSM3D
 *
 * A 3D scalar-transport blowing snow model that is spatially discretized using the finite volume method. The methods are
 * detailed in Marsh, et al (2020) but are ultimately based upon the work of Pomeroy, et al from the late 1980s and early 1990s.
 * Due to the use of the FVM, neither simplifying assumptions about wind direction nor domain rotations into the wind direction
 * are needed. This allows use where the wind ﬂow is divergent and over large extents. Erosion and deposition are computed as
 * the spatial and temporal divergence of the suspension and saltation ﬂuxes, that is, their rate of change over space and over
 * model time steps.
 *
 * The nonsteady effects of upwind fetch are represented by a downwind increase with fetch to a fully developed saturation
 * level in the saltation concentrations. This is used to calculate suspended concentrations and the increasing height of
 * the suspended snow layer with fetch. The steady‐state saltation ﬂux parameterizations (Pomeroy & Gray, 1990) are
 * used to calculate the saltation layer mass concentration based on an observed relationship between saltation trajectory
 * height and shear stress.
 *
 * Precipitation is not currently blown simultaneously
 *
 * This should be used with a wind parameterziation such as the WindNinja module.
 *
 * **Depends:**
 * - Air temp "t" [\f$ {}^\circ C \f$]
 * - Relative humidity "rh" [%]
 * - Snow water equivalent "swe" [mm]
 * - Wind speed @2m "U_2m_above_srf" [ \f$ m \cdot s^{-1} \f$ ]
 * - Wind direction @2m "vw_dir" [degrees]
 * - Wind speed @reference height "U_R" [ \f$ m \cdot s^{-1} \f$ ]
 *
 * **Provides:**
 * - Blowing snow probability. Only active if ``use_PomLi_probability:true`` "blowingsnow_probability" [0-1]
 * - Sublimation flux "Qsubl" [ \f$ kg \cdot (m^2 \cdot s)^{-1} \f$ ]
 * - Sublimation mass "Qsubl_mass" [ \f$ kg \cdot m^{-2} \f$ ]
 * - Total sublimation "sum_subl"
 * - Mass eroded or deposited during blowing snow. Positive for mass deposition, negative for mass removal. [ \f$ kg \cdot m^{-2} \f$ ]
 * - Down wind suspension flux "Qsusp" [ \f$ kg \cdot (m \cdot s)^{-1} \f$ ]
 * - Saltation flux "Qsalt" [\f$ kg \cdot (m \cdot s)^{-1} \f$]
 * - Cumulative mass erorded or desposited during model run "sum_drift"  [\f$ kg \cdot m^{-2} \f$ ]
 * - Upwind fetch if ``exp_fetch`` or ``tanh_fetch`` are used "fetch" [m]
 * - Hours since last snowfall if only ``use_PomLi_probability`` is used "p_snow_hours" [hr]
 *
 * **Parameters:**
 * - If ``enable_veg=True``, then the vegetation height "CanopyHeight" [m]
 * - If ``use_R94_lambda=True`` then leaf area index "LAI" [-]
 *
 * **Configuration:**
 *
 * \rst
 * .. code:: json
 *
 *    {
 *       "use_exp_fetch": false,
 *       "use_tanh_fetch": true,
 *       "use_PomLi_probability": false,
 *       "z0_ustar_coupling": false,
 *       "use_R94_lambda": true,
 *       "N":1,
 *       "dv":0.8,
 *       "debug_output": false,
 *       "nLayer": 10,
 *       "do_fixed_settling": false,
 *       "settling_velocity":0.5,
 *       "do_sublimation": true,
 *       "do_lateral_diff": true,
 *       "smooth_coeff":820,
 *       "min_sd_trans": 0.1,
 *       "cutoff": 0.3,
 *       "snow_diffusion_const":0.3,
 *       "rouault_diffusion_coef": false,
 *       "enable_veg": true,
 *       "iterative_subl": false,
 *
 *    }
 *
 *
 *
 * .. confval:: use_exp_fetch
 *
 *    :default: false
 *
 *    Use the exp fetch proposed by Liston Liston, G., & Sturm, M. (1998). A snow-transport model for complex terrain. Journal of Glaciology
 *
 * .. confval:: use_tanh_fetch
 *
 *    :default: true
 *
 *    Use the tanh formulation of Pomeroy & Male (1986). Has more physical basis than  ``use_exp_fetch``.
 *
 * .. confval:: use_PomLi_probability
 *
 *    :default: false
 *
 *    Considers temporal non-steady effects on blowing snow occurence via the upscaled Pomeroy and Li (2000) parameterization.
 *    Requires hours since last snowfall "p_snow_hours" [hr]. This has not been extensively tested for applicability to triangles and
 *    should be used with caution.
 *
*
 * .. confval:: z0_ustar_coupling
 *
 *    :default: false
 *
 *    To test if blowing snow occurs, a friction velocity is compared to a threshold. This friction velcoity needs a z0 estimated to be computed.
 *    If there is blowing snow, the z0 value is larger than under non-blowing snow conditions. Therefore an execution order problem occurs: should
 *    the blowing snow z0 and impact on u* be used to test for blowing snow? Setting this to true uses the blowing snow z0. Testing suggests this can
 *    overestiamte blowing snow occurence.
 *
 * .. confval:: use_subgrid_topo
 *
 *    Experimental sub-grid topographic impacts. Do not use.
 *
 * .. confval:: use_subgrid_topo_V2
 *
 *    Experimental sub-grid topographic impacts v2. Do not use.
 *
 * .. confval:: use_R94_lambda
 *
 *    :default: true
 *
 *    Estimate vegetation characteristics from LAI observation. This has been tested less extensively than the N, dv
 *    options below.   Raupach (1994) (DOI:10.1007/BF00709229)
 *

 *
 * .. confval:: N
 *
 *    :default: 1
 *
 *    If ``use_R94_lambda:false`` then use the vegetation number dnesity (number/m^2)
 *
 * .. confval:: dv
 *
 *    :default: 0.8
 *
 *    If ``use_R94_lambda:false`` then use the vegetation stalk diameter (m)
 *
 * .. confval:: debug_output
 *
 *    :default: false
 *
 *    Write extra run-time diangostics such as per-layer suspension concentration to the output file. Adds significant computtional overhead.
 *
 * .. confval:: nLayer
 *
 *    :default: 10
 *
 *    Number of vertical discretization layers. 10-15 have been found to generally work well
 *
 * .. confval:: do_fixed_settling
 *
 *    :default: false
 *
 *    If true, instead of computing a settling velocity a fixed value given by ``settling_velocity`` is used.
 *
 *
 * .. confval:: settling_velocity
 *
 *    :default: 0.5 [m/s]
 *
 *    Fixed settling velocity. Requires setting ``do_fixed_settling:true``
 *
 *
 * .. confval:: do_sublimation
 *
 *    :default: true
 *
 *    Compute sublimation losses from the saltation and suspension layers.
 *
 * .. confval:: do_lateral_diff
 *
 *    :default: true
 *
 *    Add a very small amount of lateral diffusion for numerical robustness. Should be left alone.
 *
 * .. confval:: smooth_coeff
 *
 *    :default: 820
 *
 *    Multidimensional transport equations may have spurious oscillations in the solution (Kuzmin, 2010). This term is a
 *    Laplacian smoothing (Kuzmin, 2010) term, without which oscillations between erosion and deposition may appear. This should be set
 *    to be the largest distance over which oscillations are allowed. It should be a few times the average triangle length scale (:math:`\alpha`):
 *
 *    .. math::
 *
 *       \frac{\alpha^2}{\pi^2}.
 *
 * .. confval:: min_sd_trans
 *
 *    :type: double
 *    :default: 0.1
 *
 *    Minimum snowdepth below which blowing snow saltation is not computed for this triangle. A basic form of sub-grid variability.
 *
 * .. confval:: cutoff
 *
 *    :type: double
 *    :default: 0.3
 *
 *    If (vegetation heigh - snow depth) <= cutoff, blowing snow (saltation) is inhibilited for this triangle. Used to allow
 *    some in-shrub blowing snow to occur.
 *
 * .. confval:: snow_diffusion_const
 *
 *    :type: double
 *    :default: 0.3
 *
 *    Scales the eddy diffusivity. There is some physical basis for it, but this parameter has quite a bit of
 *    variability in the literature 0.3 to 1. In PBSM3D better results have been found for values ~0.3 such that these
 *    results closely match QSusp estimates from Pomeroy, et al. This accounts for lower than reported in the literature values
 *    for the settling velocity. Thus, this parameter, for larger triangles (versus point models) accounts for spatial variability
 *    in fall velocities and turbulent diffusion. The larger triangles require setting this lower than what is generally reported
 *    for a point-scale model.
 *
 * .. confval:: rouault_diffusion_coef
 *
 *    :default: false
 *
 *    This over rides the ``snow_diffusion_const`` by using the Rouault diffusion coefficient estimation (e.g., used by Déry, et al (1999)).
 *    This approach tends to produce a diffusion coefficient that is a bit too high.
 *
 *    Rouault, M., Mestayer, P., Schiestel, R. (1991).
 *    A model of evaporating spray droplet dispersion Journal of Geophysical Research: Oceans  96(C4), 7181-7200. https://dx.doi.org/10.1029/90jc02569
 *
 * .. confval:: enable_veg
 *
 *    :default: true
 *
 *    Set to false to ignore vegetation impacts
 *
 * .. confval:: iterative_subl
 *
 *    :default: false
 *
 *    Use the Pomeroy and Li (2000) iterative solution for Schimdt's sublimation equation. This code path has not had
 *    extensive testing and should not be used at the moment.
 *
 *
 *
 * \endrst
 *
 *
 *
 * **References:**
 * - Marsh, C., Pomeroy, J., Spiteri, R., Wheater, H. (2020). A Finite Volume Blowing Snow Model for Use With Variable
 * Resolution Meshes Water Resources Research  56(2)https://dx.doi.org/10.1029/2019wr025307
 * @}
 */
class PBSM3D : public module_base
{
    REGISTER_MODULE_HPP(PBSM3D);

  public:
    PBSM3D(config_file cfg);
    ~PBSM3D();
    void run(mesh& domain);
    void init(mesh& domain);

    double nLayer;
    double susp_depth;
    double v_edge_height;

    // Beta * K, this is beta and scales the eddy diffusivity
    double snow_diffusion_const;
    double l__max;                // vertical mixing length (m)
    bool rouault_diffusion_coeff; // use the spatially variable diffusivity coefficient of Rouault 1991

    bool do_fixed_settling;   // should we have a constant settling velocity?
                              // true: constant settling velocity = settling_velocity (see below)
                              // false: use the parameterization of Pomeroy et al. (1993) and Pomeroy and Gray (1995):
                              //        In this case, the settling velocity decreases with height above the snow surface
                              //        due to a decrease with height in the mean particle size.
    double settling_velocity; // Variable used if do_fixed_settling = true
    double n_non_edge_tri;
    double eps; //lapacian smoothing epilson.
    bool do_sublimation; // should we have a sink sublimation term?
    bool do_lateral_diff; // should have lateral diffusion
    bool enable_veg;      // should we consider vegetation ?

    bool use_PomLi_probability; // Use areal Pomeroy Li 2000 probability function.
    bool use_exp_fetch;         // Enable the exp Liston 2006 fetch
    bool use_tanh_fetch;        // Enable the tanh Pomeroy and Male 1986 fetch

    bool use_subgrid_topo;    // Enable effect of subgrid topography on snow transport
    bool use_subgrid_topo_V2; // Enable effect of subgrid topography on snow transport

    bool iterative_subl; // if True, enables the iterative sublimation calculation as per Pomeroy and Li 2000
    bool use_R94_lambda; // use the ﻿Raupach 1990 lambda expression using LAI/2 instead of pomeroy stalk density



    // this is the suspension transport matrix
    double nnz; // number non zero

    // this is the drift matrix
    double nnz_drift; // number non zero

    bool debug_output;
    double cutoff; // cutoff veg-snow diff (m) that we inhibit saltation entirely

    // don't allow transport if below this threshold.
    // This gives models like snobal a chance to build up their snowpack and avoid convergence issues with thin snowcovers
    double min_sd_trans;

    // Couple the calculation of u* and z0 via the z0 value from
    // Li and Pomeroy 2000, eqn 5.
    // to modify the u* estimation instead of using a snow z0 for u* estimation
    bool z0_ustar_coupling;

    class data : public face_info
    {
      public:
        // edge unit normals
        arma::vec m[5];

        // prism areas
        double A[5];

        // face neighbors
        bool face_neigh[3];

        std::vector<double> u_z_susp; // suspension layer windspeeds
        size_t cell_local_id;

        double CanopyHeight;
        double LAI;

        // used for the pomeroy stalk formulation
        double N;  // vegetation number density
        double dv; // stalk diameter

        // saltation height
        double hs;

        bool is_edge;

        // used to flag the large vegetation areas or other via landcover types to not do any saltation at this point.
        bool saltation;

        double z0;

        double sum_drift;
        double sum_subl;
        std::vector<double> csubl; //vertical col of sublimation coeffs

        gsl_function F_fill;
        gsl_function F_fill2;
        gsl_function F_fill3;
        // gsl_function F_roots;
        // struct my_fill_topo_params params = { 1.1, 0.3, 0.6 , 0.4};
    };

    void checkpoint(mesh& domain,  netcdf& chkpt);
    void load_checkpoint(mesh& domain,  netcdf& chkpt);

private:

  // For detecting if there is suspension and/or saltation
  bool suspension_present, deposition_present;
  constexpr static double suspension_present_threshold=1e-12;
  constexpr static double deposition_present_threshold=1e-12;

  std::unique_ptr<math::LinearAlgebra::NearestNeighborProblem> deposition_NNP;
  std::unique_ptr<math::LinearAlgebra::NearestNeighborProblem> suspension_NNP;

};

/**
@}
*/
