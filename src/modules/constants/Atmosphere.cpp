#include <constants/Atmosphere.h>

namespace Atmosphere
{

// Logrithmic, following Anderson 1967? (DHSVM)  https://github.com/UW-Hydro/DHSVM/blob/master/CalcAerodynamic.c
double log_scale_wind(double u, double Z_in, double Z_out)
{
    const double Z0_SNOW  = Snow::Z0_SNOW; // Snow roughness (m)
    double Ftcr; // Multiplicative factor

    // Calc factor

    // Logarithmic wind profile assumption
    // Assuming no canopy between Z_in and Z_out

    if (Z_in > Z_out) { // If reference height is higher than measured height
        Ftcr = log((Z_out + Z0_SNOW) / Z0_SNOW) / log(Z_in / Z0_SNOW);
    } else { // If measured height is higher than reference height
        // Flip inputs, and then take 1 over Fct to get correct scaling both ways
        Ftcr = log((Z_in + Z0_SNOW) / Z0_SNOW) / log(Z_out / Z0_SNOW);
        Ftcr = 1 / Ftcr;
    }

    // Calc new Wind speed
    u = u * Ftcr;
    return u;
}

}
