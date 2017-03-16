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

#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/ell_matrix.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
// we need to define these and have them return constant values under a no-omp situation
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0

#endif

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
* Lehning_blowing_snow blowing snow impliments the blowing snow model in Alpine3D as described
 * in Lehning, et al. 2008
 * Lehning, M., Löwe, H., Ryser, M., & Raderschall, N. (2008). Inhomogeneous precipitation distribution and snow transport in steep terrain. Water Resources Research, 44(7), 1–19. http://doi.org/10.1029/2007WR006545
* Depends:
* -
*
* Provides:
* -
*/
class Lehning_blowing_snow : public module_base
{
public:
    Lehning_blowing_snow(config_file cfg);
    ~Lehning_blowing_snow();
    void run(mesh domain);
    void init(mesh domain);

    double nLayer;
    double susp_depth;
    double v_edge_height;

    double l__max; // vertical mixing length (m)
    double drift_density; // density of the transported snow to use for calculating swe depth
    bool do_vertical_advection; // should we use the discretization that includes 3D advection?

    //http://stackoverflow.com/a/4609795/410074
    template <typename T> int sng(T val) {
        return (T(0) < val) - (val < T(0));
    }

    class data : public face_info
    {
    public:
        //edge unit normals
        arma::vec m[5];

        //prism areas
        double A[5];

        //face neighbours
        bool face_neigh[3];

        //saltation height
        double hs;

        bool is_edge;

        double z0;
    };

};

/**
@}
*/
