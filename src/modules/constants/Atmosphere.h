#include <math.h>
#include <constants/Snow.h>

#pragma once

namespace Atmosphere {
    /********* Atmosphere ************/

    const double Z_U_R = 50.0; // Reference height all input wind speed forcing data are scaled to so that 1) they can be interpolated and 2) all modules can scale
    // down the reference wind speed to the height they require (i.e. through a canopy)

    const double KinVisc = 1.88e-5;     // kinematic viscosity of air (Sask. avg. value) (units ????)

    double log_scale_wind(double u, double Z_in, double Z_out);

}


