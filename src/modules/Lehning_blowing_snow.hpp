#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <cstdlib>
#include <string>
#include "math/coordinates.hpp"
#include <constants/PhysConst.h>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

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

};

/**
@}
*/
