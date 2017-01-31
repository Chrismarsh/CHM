#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include "math/coordinates.hpp"

#include <constants/PhysConst.h>
#include "constants/Atmosphere.h"

#include <cmath>


#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/bicgstab.hpp>


#include <armadillo>
#include <cstdlib>
#include <string>
//#define _USE_MATH_DEFINES
//#include <math.h>
#include <boost/math/tools/roots.hpp>
#include "math/coordinates.hpp"
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
    };

};

/**
@}
*/
