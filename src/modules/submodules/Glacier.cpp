#include "Glacier.hpp"
#include "PhysConst.h"
#include "math/Bisection.hpp"
#include <numeric>
#include <stdexcept>
#include <format>

namespace Glacier
{
    constexpr auto MM_PER_M = 1000.0;
    constexpr auto MegaGram_per_KiloGram = 1 / 1000.0;
    namespace Densification
    {
        static constexpr auto ice_density = PhysConst::rho_ice();
        
        static const double k_1(const double T)
        {
            // Empirical constants
            constexpr auto a = 11.0 * MegaGram_per_KiloGram;
            constexpr auto b = 10160.0;
            return a * std::exp(-b/(PhysConst::R() * T));

        };
        static const double Z(const double h, const double rho, const double T)
        {
            auto k = k_1(T); 
            return std::exp(ice_density * k * h) * rho / ( ice_density - rho);
        };
        static const double Z(const double h, const double rho, const double T, const double A, 
                const double rho_c)
        {
            constexpr double a = 575.0, b = 21400.0;
            auto k_2 = a * std::exp(-b/(PhysConst::R() * T)) * MegaGram_per_KiloGram;

            auto critical_depth = 1/(ice_density * k_1(T)) * 
                std::log( rho_c * ( ice_density - rho ) / (rho * (ice_density - rho_c)));

            return std::exp(ice_density * k_2 * (h - critical_depth) / std::sqrt(A)) * rho_c / (ice_density - rho_c);
        };
        static const double HerronLangwayFormula(const double Z)
        {
            return ice_density * Z / (1 + Z);
        };
        const double HerronLangway(const double depth, const double density, 
                const double glacier_temperature,const double accumulation_rate,const double critical_density)
        {
            if (density >= ice_density) [[unlikely]]
            {
                std::string err = std::format("Herron and Langway Densification formula cannot handle a density at or above {} kg/m^3. Received {} kg/m^3",ice_density,density);
                throw std::runtime_error(err);
            }

            double z;

            if (density <= critical_density)
                z = Z(depth,density,glacier_temperature);
            else
                z = Z(depth,density,glacier_temperature,accumulation_rate,critical_density);

            return HerronLangwayFormula(z);
        };

        const double Linear(const double rho, const double delta_rho)
        {
            return rho + delta_rho;
        };

    }; // Densification namespace

    static constexpr double swe_to_firn(const double h) {
        return 450.0 - 204.7 / h * (1 - std::exp(-h / 0.673));
    };

    Layer::Layer(const Units::Milimeters water_eq)
    {
        using namespace Bisection;
        auto test_function = [=,WE = water_eq](const double h) {
            auto height = std::max(h,0.6);
            auto rho = swe_to_firn(height);
            return rho * height * MM_PER_M / PhysConst::water_reference_density() - WE.value;
        };
        constexpr double low_guess = 0.6;
        constexpr double high_guess = 8.0;
		constexpr double tolerance = 1e-15;
		
        auto result = bisection(test_function,
				low_guess,high_guess,
				tolerance);

        switch (result.root) {
            case Root::Found:
                _height.value = result.value;
                _density.value = swe_to_firn(_height.value);
                break;
            case Root::DidNotConverge:
            {
                std::string err = std::format("Layer Bisection did not converge, with guess bounds ({},{}) and water equivalent {} mm",low_guess,high_guess,water_eq.value);
                throw std::runtime_error(err);
                break;
            }
            case Root::NoRootGuaranteed:
            {
                /*
                 * swe_to_firn is only well defined for snow depths above 0.6 m
                 * Therefore, for WE < rho(0.6) * 0.6, it means that the SWE value is 
                 * unsuited for this equation. Instead, take the minimal density value
                 * and recompute height
                 */
                if ( swe_to_firn(0.6) * 0.6 - water_eq.value > 0.0 )
                {
                    _density.value = swe_to_firn(0.6);
                    update_height(water_eq.value);
                }
                else
                {
                    std::string err = std::format("Layer Bisection might never converge, with guess bounds ({},{}) and water equivalent {} mm",low_guess,high_guess,water_eq.value);
                    throw std::runtime_error(err);
                }
                break;
            }
        }

    };

    Layer::Layer(const Height h, const Density rho)
    {
        _height = h;
        _density = rho;
    };

	const Layer::Height Layer::height() const { return _height;};
	const Layer::Density Layer::density() const {return _density;};

    void Layer::densifyLinear(const Density small, const Density big, const Density critical)
    {
        Layer::Density increment{0.0};
        auto WE = water_equivalent();

        if (_density.value <= critical.value)
            increment += small;
        else
            increment += big;

        _density.value = Densification::Linear(_density.value,increment.value);
        update_height(WE.value);
    };

    void Layer::densifyHerronLangway(const Units::Kelvin T, double depth,
                const double accumulation_rate, const Density critical_density)
    {
        depth += _height.value;
        auto WE = water_equivalent();

        auto new_value = Densification::HerronLangway(depth,_density.value,T.value,accumulation_rate,critical_density.value);
        if (new_value > _density.value)
        {
            _density.value = new_value;
            update_height(WE.value);
        }
    };

    void Layer::update_height(const double old)
    {
        
        _height.value = old  * PhysConst::water_reference_density() / (_density.value * MM_PER_M);
    };

    void Layer::remove(const Units::Milimeters amount)
    {
        auto current = water_equivalent();
        if (amount >= current)
            throw std::invalid_argument("Layer::remove amount must not be in excess of or equal to current");

        // Not the same as update_height, this is a flat height reduction without changing density
        current -= amount;
        _height.value = current.value * PhysConst::water_reference_density() / 
            ( _density.value * MM_PER_M );
    };

    const Units::Milimeters Layer::water_equivalent() const
    {
        return Units::Milimeters{_height.value * _density.value / PhysConst::water_reference_density() * MM_PER_M};
    };

    LayeredFirn::LayeredFirn(const Params* _p) : p(_p) {};

    LayeredFirn::LayeredFirn(const Params* _p,const std::vector<Layer>& L) : p(_p) {
        layers.assign(L.begin(),L.end());
    };

    LayeredFirn::LayeredFirn(const Params* _p, const std::deque<Layer>& L) : p(_p), layers(L) {};

    const Units::Milimeters LayeredFirn::water_equivalent() const
    {
        if (layers.empty())
            return Units::Milimeters{0.0};
        
        auto sum = std::accumulate(layers.begin(),layers.end(),0.0,
                [](double sum, const Layer& layer) {
                    return sum + layer.water_equivalent().value;
                    });

        return Units::Milimeters{sum};
    };

    void LayeredFirn::accumulate(const Units::Milimeters WE)
    {
        Layer new_layer(WE);
        layers.push_back(new_layer);
    };

    const MeltInfo LayeredFirn::melt(const Units::Milimeters input)
    {
		Units::Milimeters melt{0.0};
		Units::Milimeters init_layer_WE{0.0};
        while (input.value > 0.0 && !layers.empty())
        {
            Layer& current = layers.back();

            if ( input.value < current.water_equivalent().value )
            {
				init_layer_WE = current.water_equivalent();
                current.remove(input);
                melt += input;
                input.value = 0.0;
				if (current.water_equivalent().value < 1e-12)
					layers.pop_back();
            }
            else
            {
                const auto value = current.water_equivalent();
                melt += value;
                input -= value;
                layers.pop_back();
            };
        };
		
		if (input.value > 0.0 && input.value < 1e-12)
		{
			input.value = 0.0;
			layers.pop_back();
		};

        return MeltInfo {
            melt,
            input			
        };

    };

    const LayeredFirn::LayerContainer& LayeredFirn::get_layers() const
    {
        return layers;
    };

	std::optional<Units::Milimeters> LayeredFirn::convert_to_ice()
    {
        std::optional<Units::Milimeters> result;

        if (layers.empty())
            return result;

        if (layers.front().density().value > p->firn_to_ice_density)
        {
            result.emplace(layers.front().water_equivalent());
            layers.pop_front();
        }
        return result;
    };

    Ice::Ice(const Params* _p) : p(_p) {};

    Ice::Ice(const Params* _p, const Units::Milimeters WE) : p(_p), _water_equivalent(WE) {};

    void Ice::accumulate(const Units::Milimeters WE)
    {
        _water_equivalent += WE;
    };

    const MeltInfo Ice::melt(const Units::Milimeters input)
    {
        auto melted = Units::Milimeters{0.0};

        if (input > _water_equivalent)
        {
            melted += _water_equivalent;
			_water_equivalent.value = 0.0;
        }
        else
		{
            melted += input;
			_water_equivalent -= input;
		}

        Units::Milimeters remaining_energy = input - melted;
        remaining_energy.value = std::max(remaining_energy.value,0.0);

        return MeltInfo{melted,
                remaining_energy};
    };

    const Units::Milimeters Ice::water_equivalent() const
    {
        return _water_equivalent;
    }

	namespace Melt
	{
        Units::Milimeters convert_to_mass(const Units::Watts_per_m2 energy,const size_t dt)
        {
            constexpr auto thermal_factor = 0.95;
            return energy.value * dt / (PhysConst::Lf() * PhysConst::water_reference_density() * thermal_factor);
        };

		MeltScenario get_scenario(const State& s,
				Units::Milimeters melt_energy,const Units::Milimeters swe)
		{
			if (swe.value > 0.0 || melt_energy.value == 0.0)
				return MeltScenario::NoMelt;

			if (s.firn.water_equivalent().value > 0.0
					&& s.ice.water_equivalent().value == 0.0)
				return MeltScenario::FirnMelt;

			if (s.firn.water_equivalent().value == 0.0
					&& s.ice.water_equivalent().value > 0.0)
				return MeltScenario::IceMelt;

			if (s.firn.water_equivalent().value > 0.0
					&& s.ice.water_equivalent().value > 0.0)
			{
				if ( melt_energy.value > s.firn.water_equivalent().value)
					return MeltScenario::FirnAndIceMelt;
				else
					return MeltScenario::FirnMelt;
			}

			if (s.total_water_equiv().value == 0.0)
				return MeltScenario::NoMelt;

			std::string err = 
				std::format("Unhandled melt scenario in glacier submodule for SWE of {} mm, firn of {} mm, and ice of {} mm",swe.value,s.firn.water_equivalent().value,s.ice.water_equivalent().value);
			throw std::logic_error(err);
		};

		std::optional<MeltInfo> compute_melt(MeltScenario scenario,State& s,const Units::Milimeters melt_mass_max)
		{
			std::optional<MeltInfo> info;
			switch(scenario)
			{
				case MeltScenario::NoMelt:
					// No melt, so info remains uninitialized
					break;
				case MeltScenario::FirnMelt:
					info.emplace(s.firn.melt(melt_mass_max));
					break;
				case MeltScenario::FirnAndIceMelt: {
					info.emplace(s.firn.melt(melt_mass_max));
					MeltInfo ice_info = s.ice.melt(info->remaining_energy);
					info->melt += ice_info.melt;
					info->remaining_energy += ice_info.remaining_energy;
				}
					break;
				case MeltScenario::IceMelt:
					info.emplace(s.ice.melt(melt_mass_max));
					
			}
			return info;
					
		};


	};

    namespace Updates
    {
        void swe_to_firn(const Params& p, LayeredFirn& firn, Units::Milimeters& swe)
        {
            if (swe.value == 0.0)
                return;

            firn.accumulate(swe);

            swe.value = 0.0;  
             
        };

        void firn_to_ice(const Params& p, LayeredFirn& firn, Ice& ice)
        {
            auto new_ice = firn.convert_to_ice();

            if (new_ice)
                ice.accumulate(*new_ice);

        };
    }

	State::State(const Params* p) : firn(p), ice(p) {};

	State::State(const Params* p, const Ice& i) : firn(p), ice(i) {};

	State::State(const LayeredFirn& f,const Ice& i) : firn(f), ice(i) {};

	const Units::Milimeters State::total_water_equiv() const
	{
		auto result = ice.water_equivalent();
		result += firn.water_equivalent();

		return result;
	};

	const Units::Metres State::total_depth() const
	{
		using namespace PhysConst;
		Units::Metres sum{water_reference_density() / rho_ice() *
			ice.water_equivalent().value / MM_PER_M};

		auto& layers = firn.get_layers();

		for (auto layer : layers)
		{
			sum += layer.height();
		}

		return sum;
	};
};
