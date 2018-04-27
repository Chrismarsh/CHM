/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

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
    const double rho_ice = 917;     //{density of ice (kg/m3)} 


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

}


