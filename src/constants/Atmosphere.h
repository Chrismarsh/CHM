#include <math.h>
#include <constants/Snow.h>

#pragma once

namespace Atmosphere {
    /********* Atmosphere ************/

    const double Z_U_R = 50.0; // Reference height all input wind speed forcing data are scaled to so that 1) they can be interpolated and 2) all modules can scale
    // down the reference wind speed to the height they require (i.e. through a canopy)

    const double KinVisc = 1.88e-5; // kinematic viscosity of air (Sask. avg. value) (units ????)

    double log_scale_wind(double u, double Z_in, double Z_out, double snowdepthavg, double z0=Snow::Z0_SNOW);

    // Inoue E (1963) On the turbulent structure of air flow within crop canopies. J Meteorol Soc Jpn 41:317–326
    double exp_scale_wind(double u, double Z_in, double Z_out, const double alpha);

}


