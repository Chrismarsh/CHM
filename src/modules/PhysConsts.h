
#ifndef __PHYSCONSTS_H__
#define __PHYSCONSTS_H__

namespace PhysConsts {

    //TODO: CHeck values against ones in CRHM (waiting for values from Tom Brown)

    //const double Tm = mio::Cst::t_water_freezing_pt; // 273.15; // (K)
    const double emiss_c    = 0.5; // Emissitivity of canopy TODO: GUESS VALUE: Need to get from lookup table depending on plant type
    const double sbc        = mio::Cst::stefan_boltzmann;
    const double Rgas       = mio::Cst::gaz_constant_dry_air; // (W m-2 K-4) TODO: Check this is the right constant
    const double kappa      = 0.4; // Von Karman constant
    const double Ls         = mio::Cst::l_water_sublimation; // (J kg-1) TODO: What is this?
    const double emiss      = 0.99; // Long-wave Emissivity of snow at nadir
    const double Cp         = mio::Cst::specific_heat_air; // 1004.67; // (J K-1), see Stull "Meteorology for scientists and engineers" p44 TODO: CHECK CRHM value

    // From PBSM Pomeroy et al. 1993
    const double M          = 18.01; // M is the molecular weight of water (18.01 kg kmol-l); TODO: CHECK CRHM value
    const double R          = mio::Cst::gaz_constant; // R is the universal gas constant (8313 J mol-1 K-1); FYI: unit error in Pomeroy et al. 1993 TODO: CHECK CRHM value
    const double DICE       = nullptr; // TODO: What is this?


}