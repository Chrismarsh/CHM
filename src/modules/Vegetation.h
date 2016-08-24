
#pragma once

namespace Vegetation {
    /********* Vegetation ************/
    const double emiss_c    = 0.96;     // emissivity of a canopy
    const double alb_c      = 0.1;      // "canopy albedo" 0.05-0.2 TODO: take from vege look up table
    const double ks         = 0.0114;   // snow shape coefficient for jack pine TODO: take from vege look up table
    const double Fract      = 0.37;     // fractal dimension of intercepted snow TODO: take from vege look up table
}

#endif //CHMPRJ_VEGETATION_H
