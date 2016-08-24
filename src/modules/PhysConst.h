#include <meteoio/MeteoIO.h>

#pragma once

namespace PhysConst {
    /********* Physical Constants ************/
    const double sbc = mio::Cst::stefan_boltzmann; // 5.670373e-8 (W m-2 K-4)
    const double Rgas = mio::Cst::gaz_constant_dry_air; // (W m-2 K-4)
    const double kappa = 0.4; // Von Karman constant
    const double Ls = mio::Cst::l_water_sublimation; // 2.838e6; // Latent heat of sublimation // ((J kg-1)
    const double Cp = mio::Cst::specific_heat_air; // 1004.67; // (J K-1), see Stull "Meteorology for scientists and engineers" p44
    const double Ci = mio::Cst::specific_heat_ice; //  = 2100.0 (J K-1), at 0C

    // From PBSM Pomeroy et al. 1993
    const double M = 18.01; // M is the molecular weight of water (18.01 kg kmol-l);
    const double R = mio::Cst::gaz_constant; // R is the universal gas constant (8313 J mol-1 K-1); FYI: unit error in Pomeroy et al. 1993
    const double DICE = 900;     //{density of ice, assumed equal to blowing snow part (kg/m3)} TODO: This is incorrect (ask John why not 917)


    // From CRHM_constants in common.h
    const double Cs = 1.28E+06; // (J/m3/K) volumetric heat capacity of soil
    //const double Cp = 1005;     // (J/kg/K) volumetric heat capacity of dry air
    //const double Rgas = 287.0;  // Gas constant for dry air (J/kg/K)
    //const double Tm = 273.15;   // Melting point (K)

    //const double Ls = 2.845e6;  // Latent heat of sublimation (J/kg)
    const double Lv = 2.50e6;  // Latent heat of vaporization (J/kg)
    const double Lf = 0.334e6;  // Latent heat of fusion (J/kg)
    //const double kappa = 0.4;

    //const double sbc = 5.67E-8; // Stephan-Boltzmann constant W/m^2/k4
    const double SB = 4.899e-09; // Stephan-Boltzmann constant MJ/m^2-d

    const double em = 0.622;     //




    /*
     * namespace PBSM_constants {

  const float rho = 1.23;     // (kg/m3) density of dry air
  const float Qstar = 120;    //{Solar Radiation Input}
  const float M = 18.01;      //{molecular weight of water (kg/kmole)}
  const float R = 8313;       //{universal gas constant (J/(kmole K))}
  const float LATH = 2.838E6; //{latent heat of sublimation (J/kg) List 1949}
  const float DICE = 900;     //{density of ice, assumed equal to blowing snow part (kg/m3)}
  const float ZD = 0.3;       //{height of boundary-layer at xd (m) Takeuchi (1980)}
  const float XD = 300;       //{Fetch to develop b-l to ZD (m)}
  const float g = 9.80;       //{m/s2}
  const float Beta = 170;     // Beta ratio of element to surface drag for vegetation Cr/Cs
  const float C1 = 2.8;       //{2.3}
  const float C2 = 1.6;
  const float C3 = 4.2;       //{3.25} {e = 1/(C3*Ustar)}
  const float KARMAN = 0.4;
  const float KARMAN2 = 0.16;
     */

}


#endif //CHMPRJ_PHYSCONST_H
