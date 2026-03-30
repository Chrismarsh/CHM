#pragma once

#include "Atmosphere.h"
#include "base_step.hpp"
#include "PhysConst.h"
#include <concepts>

namespace Units = PhysConst::units;

namespace katabatic_melt_energy
{
    struct Params
    {
        const double heat_capacity_air = PhysConst::Cp(); 
        const double molecular_wt_ratio = PhysConst::em();
        const double g = PhysConst::g();  

        double prandtl = 5.0;
        double k = 4e-4;
        double k2 = 1.0;
		size_t seconds_per_step = 3600.0;
    }; 

    const Units::m_per_s bulk_coefficient( const Params&, const Units::TempDiff deficit,
            const Units::LapseRateSI gamma,
            const Units::Kelvin glacier_temperature);

    const Units::Watts_per_m2 sensible_heat(const Params&, const Units::m_per_s coeff,const Units::TempDiff deficit,const Units::DensitySI air_density);

    const Units::Watts_per_m2 latent_heat(const Params&, const Units::m_per_s coeff,const Units::Pa vapour_pressure_deficit,
            const Units::Kelvin glacier_temperature, const Units::DensitySI air_density, const Units::Pa air_pressure);
	
	template<typename T>
	concept KatabaticData = requires(T& t)
    {
		{ t.glacier_temperature() } -> std::same_as<Units::Kelvin>;
		{ t.air_temperature() } -> std::same_as<Units::Kelvin>;
		{ t.air_pressure() } -> std::same_as<Units::Pa>;
		{ t.vapour_pressure() } -> std::same_as<Units::Pa>;
		{ t.vapour_pressure_surface() } -> std::same_as<Units::Pa>;
		{ t.lapse_rate() } -> std::same_as<Units::LapseRateSI>;

		{ t.latent_heat(std::declval<double>()) } -> std::same_as<void>;
		{ t.sensible_heat(std::declval<double>()) } -> std::same_as<void>;
	};

	template<KatabaticData Data>
	class Model : public base_step<Model<Data>,Data>
	{
		Params p;
	public:
		void execute_impl(Data& d) const
		{
			auto air_pressure = d.air_pressure();
			auto vapour_pressure = d.vapour_pressure();
			auto glacier_temperature = d.glacier_temperature();
			auto air_temperature = d.air_temperature();
			auto lapse_rate = d.lapse_rate();

			auto rho_air = Atmosphere::air_density(air_pressure,glacier_temperature,vapour_pressure);

			Units::TempDiff T_deficit{air_temperature.value - glacier_temperature.value};

			auto K = bulk_coefficient(p,T_deficit,lapse_rate,glacier_temperature);

			auto Q_sensible = sensible_heat(p,K,T_deficit,rho_air);

			Units::Pa vapour_pressure_deficit{vapour_pressure.value - d.vapour_pressure_surface().value};

			auto Q_latent = latent_heat(p,K,vapour_pressure_deficit,glacier_temperature,rho_air,air_pressure);

			d.latent_heat(Q_latent.value);
			d.sensible_heat(Q_sensible.value);
		};

		Params& get_params() { return p; };
	};
};
