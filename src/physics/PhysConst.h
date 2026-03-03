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
#include <type_traits>


namespace PhysConst {
	namespace units
	{
		template<typename T,typename Derived>
		struct base { 
			T value; 
			
			static_assert(std::is_arithmetic_v<T>,"PhysConst::units::base template must use a built-in integral or floating-point type");

                // Prevent derived classes from getting their own synthesized comparisons

			void operator+=(const Derived& other) 
			{
				static_cast<Derived&>(*this).value += other.value;
			}
			void operator-=(const Derived& other)
			{
				static_cast<Derived&>(*this).value -= other.value;
			}
            friend auto operator<=>(const Derived& lhs,const Derived& rhs)
            {
                return lhs.value <=> rhs.value; 
            } 
            friend bool operator==(const Derived& lhs,const Derived& rhs)
            {
                return lhs.value == rhs.value; 
            }
		};

        // Temperature
		static inline constexpr double freeze_point = 273.15;
		struct Celsius;
		struct Kelvin 
			: public base<double,Kelvin> {
				Kelvin(const Celsius& c);
				explicit Kelvin(const double in) { value = in; };
			};
		struct Celsius
			: public base<double,Celsius> {
				Celsius(const Kelvin& k)
				{
					value = k.value - freeze_point;
				};
				explicit Celsius(const double in) {value = in;};
			};
		inline Kelvin::Kelvin(const Celsius& c) { value = c.value + freeze_point; };
		struct TempDiff
			: public base<double,TempDiff> {
				TempDiff(const Kelvin &c) {value = c.value; };
				TempDiff(const Celsius &c) {value = c.value; };
				explicit TempDiff(const double in) { value = in; };
			};


        //depths
        struct Milimeters
            : public base<double,Milimeters> {};
        struct Metres
            : public base<double,Metres> {};

        // mass
        struct Kg_per_m3
            : public base<double,Kg_per_m3> {};
        using DensitySI = Kg_per_m3;


        // pressure
        struct Pa
            : public base<double,Pa> {};
        
        // Lapse rate
        struct Degree_per_metre
            : public base<double,Degree_per_metre> {};
        using LapseRateSI = Degree_per_metre;

		// Energy rate
		struct Watts_per_m2
			: public base<double,Watts_per_m2> {};
	};
    
     
    /********* Physical Constants ************/
    // Stefan-Boltzmann constant (W m-2 K-4)
	inline const double sbc() { return mio::Cst::stefan_boltzmann; }
    
    // Gas constant for dry air (J kg-1 K-1)
    inline const double RgasDry() { return mio::Cst::gaz_constant_dry_air; }
	// Gas constant for water vapour (J kg-1 K-1) 
    inline const double RgasVapour() { return mio::Cst::gaz_constant_water_vapor; }
    // Von Karman constant (dimensionless)
    constexpr double kappa() { return 0.4; }
    
    // Latent heat of sublimation (J kg-1)
    inline const double Ls() { return mio::Cst::l_water_sublimation; }
    
    // Specific heat of dry air at constant pressure (J kg-1 K-1)
    inline const double Cp() { return mio::Cst::specific_heat_air; }
    
    // Specific heat of ice (J kg-1 K-1)
    inline const double Ci() { return mio::Cst::specific_heat_ice; }
    
    // Molecular weight of water (kg kmol-1)
    constexpr double M() { return 18.01; }
    
    // Universal gas constant (J mol-1 K-1)
    inline const double R() { return mio::Cst::gaz_constant; }
    
    // Density of ice (kg m-3)
    constexpr double rho_ice() { return 917.0; }
    
    // Volumetric heat capacity of soil (J m-3 K-1)
    constexpr double Cs() { return 1.28E+06; }
    
    // Latent heat of vaporization, constant value (J kg-1)
    constexpr double Lv() { return 2.501e6; }
    
    // Latent heat of vaporization with temperature dependence (J kg-1)
    const double Lv(const units::Celsius T);
    
    // Latent heat of fusion (J kg-1)
    constexpr double Lf() { return 0.334e6; }
    
    // Stefan-Boltzmann constant, day units (MJ m-2 d-1)
    constexpr double SB() { return 4.899e-09; }
    
    // Ratio of molecular weights of water vapor to dry air (dimensionless)
    constexpr double em() { return 0.622; }

    // Acceleration due to gravity (m s-2)
    constexpr auto g() {return 9.81;}
    
    // Reference density of water (kg m-3)
    constexpr double water_reference_density() { return 1000.0; }

}
