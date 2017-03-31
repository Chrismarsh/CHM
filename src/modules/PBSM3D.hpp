#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "interpolation.h"

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

#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/linalg/lu.hpp>
#include <viennacl/ell_matrix.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>

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
    bool do_vertical_advection; // should we use the discretization that includes 3D advection?
    double settling_velocity;
    double n_non_edge_tri;
    double eps; //lapacian smoothing epilson.
    bool limit_mass; // do not saltate more snow than what exists in a cell.
    bool do_sublimation; // should we

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

        size_t cell_id;


        //saltation height
        double hs;

        bool is_edge;

        double z0;
    };

};

/**
@}
*/
