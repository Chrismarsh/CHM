
#ifndef __PHYSCONSTS_H__
#define __PHYSCONSTS_H__

namespace PhysConsts {

    //TODO: CHeck values against ones in CRHM (waiting for values from Tom Brown)

    //const double Tm = mio::Cst::t_water_freezing_pt; // 273.15; // (K)
    const double sbc        = mio::Cst::stefan_boltzmann; // 5.670373e-8 (W m-2 K-4)
    const double Rgas       = mio::Cst::gaz_constant_dry_air; // (W m-2 K-4)
    const double kappa      = 0.4; // Von Karman constant
    const double Ls         = mio::Cst::l_water_sublimation; // 2.838e6; // ((J kg-1)
    //const double emiss      = 0.99; // Long-wave Emissivity of snow at nadir
    const double Cp         = mio::Cst::specific_heat_air; // 1004.67; // (J K-1), see Stull "Meteorology for scientists and engineers" p44 TODO: CHECK CRHM value

    // From PBSM Pomeroy et al. 1993
    const double M          = 18.01; // M is the molecular weight of water (18.01 kg kmol-l); TODO: CHECK CRHM value
    const double R          = mio::Cst::gaz_constant; // R is the universal gas constant (8313 J mol-1 K-1); FYI: unit error in Pomeroy et al. 1993 TODO: CHECK CRHM value
    const double DICE       = 0.5; // TODO: What is this?


    // From CRHM_constants in common.h
    const double Cs = 1.28E+06; // (J/m3/K) volumetric heat capacity of soil
    //const double Cp = 1005;     // (J/kg/K) volumetric heat capacity of dry air
    //const double Rgas = 287.0;  // Gas constant for dry air (J/kg/K)
    //const double Tm = 273.15;   // Melting point (K)

    //const double Ls = 2.845e6;  // Latent heat of sublimation (J/kg)
    const double Lv = 2.50e6 ;  // Latent heat of vaporization (J/kg)
    const double Lf = 0.334e6;  // Latent heat of fusion (J/kg)
    //const double kappa = 0.4;

    //const double sbc = 5.67E-8; // Stephan-Boltzmann constant W/m^2/k4
    const double SB = 4.899e-09; // Stephan-Boltzmann constant MJ/m^2-d

    const double emiss = 0.985; // emissivity of the atmosphere and snowpack
    const double emiss_c = 0.96; // emissivity of the canopy
    const double em = 0.622;     //
    
}

#endif