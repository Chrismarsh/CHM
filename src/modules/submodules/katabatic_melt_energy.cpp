#include "katabatic_melt_energy.hpp"
#include "PhysConst.h"
#include <format>
#include <stdexcept>

namespace katabatic_melt_energy
{
    const water_flux<FluxType::Default> bulk_coefficient(const Params& p, const Units::TempDiff deficit,
            const Units::LapseRateSI gamma,
            const Units::Kelvin glacier_temperature)
    {
		if (gamma == Units::LapseRateSI{0.0} ||
				p.prandtl == 0.0 ||
				glacier_temperature == Units::Kelvin{0.0})
		{
			std::string err = std::format("Divide by zero in katabatic_melt_energy::bulk_coefficient for gamma = {} K/m, Pr = {}, T = {} K",gamma.value,p.prandtl,glacier_temperature.value);
			throw std::runtime_error(err);
		}
		const auto result = water_flux<FluxType::Default>::from_m_per_s(p.k * 
				std::pow(p.k2,2.0) * deficit.value * 
				std::sqrt(p.g / 
					(glacier_temperature.value * gamma.value * p.prandtl)));
		return result;
    };

    const water_flux<FluxType::latent> sensible_heat(const Params& p, const water_flux<FluxType::Default> coeff,const Units::TempDiff deficit,const Units::DensitySI air_density)
    {
        auto m = water_flux<FluxType::latent>::from_W_per_m_squared(air_density.value * p.heat_capacity_air 
                * coeff.m_per_s() * deficit.value);

        return m;
    };

    const water_flux<FluxType::latent> latent_heat(const Params& p, const water_flux<FluxType::Default> coeff,const Units::Pa vapour_pressure_deficit,
            const Units::Kelvin glacier_temperature, const Units::DensitySI air_density)
    {
		Units::Celsius Temp{glacier_temperature.value - 273.15};
        auto m = water_flux<FluxType::latent>::from_W_per_m_squared(p.molecular_wt_ratio * air_density.value * 
                PhysConst::Lv(Temp) * coeff.m_per_s() * vapour_pressure_deficit.value);

        return m;

    };

};
