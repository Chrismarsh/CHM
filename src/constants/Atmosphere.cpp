#include <constants/Atmosphere.h>

namespace Atmosphere
{

    // Logrithmic, assuming no snow cover and no canopy bewteen Z_in and Z_out
    // Because filter is before runtime, we do not know what the snowdepth will be, thus this
    // introduces some error latter when wind is scaled down taking into account the snowdepth.
    double log_scale_wind(double u, double Z_in, double Z_out, double snowdepthavg, double z0)
    {
        double Z0_SNOW  = z0; //Snow::Z0_SNOW; // Snow roughness (m)

        u = u * log((Z_out - (snowdepthavg + Z0_SNOW)) / Z0_SNOW) / log((Z_in - (snowdepthavg + Z0_SNOW)) / Z0_SNOW);
        return u;
    }

    // Exponential following Inoue (1963)
    double exp_scale_wind(double u, double Z_in, double Z_out, const double alpha)
    {

        u = u * exp(alpha*(Z_out/Z_in-1));
        return u;
    }

}
