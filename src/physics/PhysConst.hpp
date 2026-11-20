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
#pragma once

#include <meteoio/MeteoIO.h>

namespace PhysConst {
    /********* Physical Constants ************/
    // Stefan-Boltzmann constant (W m-2 K-4)
    const double sbc = mio::Cst::stefan_boltzmann;
    
    // Gas constant for dry air (J kg-1 K-1)
    const double RgasDry = mio::Cst::gaz_constant_dry_air;
	// Gas constant for water vapour (J kg-1 K-1) 
    const double RgasVapour = mio::Cst::gaz_constant_water_vapor;
    // Von Karman constant (dimensionless)
    constexpr double kappa = 0.4;
    
    // Latent heat of sublimation (J kg-1)
    const double Ls = mio::Cst::l_water_sublimation;
    
    // Specific heat of dry air at constant pressure (J kg-1 K-1)
    const double Cp = mio::Cst::specific_heat_air;
    
    // Specific heat of ice (J kg-1 K-1)
    const double Ci = mio::Cst::specific_heat_ice;
    
    // Molecular weight of water (kg kmol-1)
    constexpr double M = 18.01;
    
    // Universal gas constant (J mol-1 K-1)
    const double R = mio::Cst::gaz_constant;
    
    // Density of ice (kg m-3)
    constexpr double rho_ice = 917.0;

    // Volumetric heat capacity of soil (J m-3 K-1)
    constexpr double Cs = 1.28E+06;
    
    // Latent heat of vaporization, constant value (J kg-1)
    constexpr double Lv =  2.501e6;
    
    // Latent heat of fusion (J kg-1)
    constexpr double Lf = 0.334e6;
    
    // Stefan-Boltzmann constant, day units (MJ m-2 d-1)
    constexpr double SB = 4.899e-09;
    
    // Ratio of molecular weights of water vapor to dry air (dimensionless)
    constexpr double em = 0.622;

    // Acceleration due to gravity (m s-2)
    constexpr auto g = 9.81;
    
    // Reference density of water (kg m-3)
    constexpr double water_reference_density = 1000.0;

}
