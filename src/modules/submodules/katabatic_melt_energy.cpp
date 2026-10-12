#include "katabatic_melt_energy.hpp"
#include "PhysConst.h"
#include <format>
#include <stdexcept>

namespace katabatic_melt_energy
{
    const Units::m_per_s bulk_coefficient(const Params& p, const Units::TempDiff deficit,
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
		auto result = Units::m_per_s{p.k * 
				std::pow(p.k2,2.0) * deficit.value * 
				std::sqrt(p.g / 
					(glacier_temperature.value * gamma.value * p.prandtl))};
		//result.value = (0.01 + result.value)/2;
		return result;
    };

    const Units::Watts_per_m2 sensible_heat(const Params& p, const Units::m_per_s coeff,const Units::TempDiff deficit,const Units::DensitySI air_density)
    {
        auto m = Units::Watts_per_m2{air_density.value * p.heat_capacity_air 
                * coeff.value * deficit.value};

        return m;
    };

    const Units::Watts_per_m2 latent_heat(const Params& p, const Units::m_per_s coeff,const Units::Pa vapour_pressure_deficit,
            const Units::Kelvin glacier_temperature, const Units::DensitySI air_density,const Units::Pa air_pressure)
    {
		Units::Celsius Temp{glacier_temperature.value - 273.15};
        auto m = Units::Watts_per_m2{p.molecular_wt_ratio * air_density.value * 
                PhysConst::Lv(Temp) * coeff.value * vapour_pressure_deficit.value / air_pressure.value};

        return m;

    };

};
