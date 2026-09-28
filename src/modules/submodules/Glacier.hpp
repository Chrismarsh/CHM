#pragma once

#include "base_step.hpp"
#include <concepts>
#include <deque>
#include <optional>
#include "PhysConst.h"

namespace Units = PhysConst::units;
namespace Glacier
{
    enum class DensifyVersion {
        HerronLangway,
        Linear
    };

	struct Params
	{
		const double water_density =
			PhysConst::water_reference_density();
		const double heat_capacity_air = 
			PhysConst::Cp();
        const double latent_heat_vapour =
            PhysConst::Lv();
		const double gas_constant_dry =
			PhysConst::RgasDry();
		const double gas_constant_vapour =
			PhysConst::RgasVapour();
        DensifyVersion densify_version = DensifyVersion::HerronLangway;
        double small_increment{25.0};
        double big_increment{50.0};
        double critical_density{550.0};
		double firn_to_ice_density{830.0};
		double thermal_factor{0.95};
        size_t seconds_per_step{3600};
	};

    namespace Densification
    {
        // Low density
        const double HerronLangway(const double depth, const double density, 
                const double glacier_temperature);
        // Above density
        const double HerronLangway(const double depth, const double density,
                const double glacier_temperature, const double accumulation_rate, 
                const double critical_density = 550.0);

        const double Linear(const double density,const double density_increment);
    };


    class Layer
    {
    public:
        struct Height : public Units::Metres {};
        struct Density : public Units::DensitySI {}; 

	    Layer(const Units::Milimetres water_equivalent);
        Layer(const Height,const Density);

        const Units::Milimetres water_equivalent() const;
        void densifyLinear(const Density small, const Density big,const Density critical);
        void densifyHerronLangway(const Units::Kelvin, double depth,
                const double accumulation_rate, const Density critical_density);
        void remove(const Units::Milimetres);
	
		const Height height() const;
		const Density density() const;
    
    private:
        Height _height;
        Density _density;

        void update_height(const double old);
        
    };

    struct MeltInfo
    {
        Units::Milimetres melt;
        Units::Milimetres remaining_energy;
    };

    class LayeredFirn
    {
    public: 
        using LayerContainer = std::deque<Layer>;

    private:
        const Params* const p;
        LayerContainer layers;

    public:
        explicit LayeredFirn(const Params*);
        LayeredFirn(const Params*,const std::vector<Layer>&);
        LayeredFirn(const Params*,const LayerContainer&);
        std::optional<Units::Milimetres> convert_to_ice();

        void accumulate(const Units::Milimetres);
        
        const MeltInfo melt(const Units::Milimetres);
        const Units::Milimetres water_equivalent() const;
        const LayerContainer& get_layers() const;
    };
    
    class Ice
    {
        const Params* const p;
        Units::Milimetres _water_equivalent{0.0};
    public:
        explicit Ice(const Params*);
        Ice(const Params*, const Units::Milimetres);

        void accumulate(const Units::Milimetres);
        const MeltInfo melt(const Units::Milimetres);
        const Units::Milimetres water_equivalent() const;
    };

	struct State
    {
        LayeredFirn firn;
        Ice ice;

		State() = delete;
		State(const Params* p);
		State(const Params* p, const Ice& i);
		State(const Params* p, const LayeredFirn& f);
		State(const LayeredFirn& f, const Ice& i);
		const Units::Milimetres total_water_equiv() const;
		const Units::Metres total_depth() const;
    };

	namespace Melt
	{
		enum class MeltScenario {
			NoMelt, FirnMelt, FirnAndIceMelt, IceMelt};

        Units::Milimetres convert_to_mass(const Units::Watts_per_m2,const Params&);

		MeltScenario get_scenario(const State&,const Units::Milimetres melt_energy,const Units::Milimetres swe);

		std::optional<MeltInfo> compute_melt(MeltScenario,State&,const Units::Milimetres);
	};

    namespace Updates
    {
        void swe_to_firn(const Params&, LayeredFirn&, Units::Milimetres& swe);
        void firn_to_ice(const Params&, LayeredFirn&, Ice&);
    };

    template<class T>
    concept GlacierData = requires(T& t)
    {
		{ t.get_state() } -> std::same_as<State&>;
        { t.swe() } -> std::same_as<const Units::Milimetres>;
        { t.melt_energy() } -> std::same_as<const Units::Watts_per_m2>;
		{ t.glacier_temperature() } -> std::same_as<const Units::Celsius>;
        { t.update_now() } -> std::same_as<bool>;
				
    	{ t.glacier_water_equivalent(std::declval<const double>()) } -> std::same_as<void>;
        { t.total_depth(std::declval<const double>()) } -> std::same_as<void>;
        { t.firn_melt(std::declval<const double>()) } -> std::same_as<void>;
        { t.ice_melt(std::declval<const double>()) } -> std::same_as<void>;
    };

    template<GlacierData data>
    class Model : public base_step<Model<data>,data>
    {
    public:
        void execute_impl(data& d) const
        {
            auto& s = d.get_state();
			auto SWE = d.swe();
            Units::Milimetres melt_energy_mm_per_dt = Melt::convert_to_mass(d.melt_energy(),p);

            if (melt_energy_mm_per_dt.value > 0.0)
            {
                using namespace Melt;
                auto scenario = get_scenario(s,
                        melt_energy_mm_per_dt,SWE);

                std::optional<State> state_copy;

                if (scenario == MeltScenario::FirnAndIceMelt)
                    state_copy.emplace(s);

                auto info = compute_melt(scenario,s,melt_energy_mm_per_dt);

                switch (scenario) {
                    case MeltScenario::NoMelt:
                        d.firn_melt(0.0);
                        d.ice_melt(0.0);
                        break;
                    case MeltScenario::FirnMelt:
                        d.firn_melt(info->melt.value);
                        d.ice_melt(0.0);
                        break;
                    case MeltScenario::IceMelt:
                        d.ice_melt(info->melt.value);
                        d.firn_melt(0.0);
                        break;
                    case MeltScenario::FirnAndIceMelt:
                        d.firn_melt(state_copy->firn.water_equivalent().value);

                        d.ice_melt(state_copy->ice.water_equivalent().value
                                - s.ice.water_equivalent().value);
                        break;
                }
            }

            if (d.update_now())
            {
				using namespace Updates;

				swe_to_firn(p,s.firn,SWE);

				firn_to_ice(p,s.firn,s.ice);
            }

			d.glacier_water_equivalent(s.total_water_equiv().value);
			d.total_depth(s.total_depth().value);
        };
        Params& get_params()
		{
			return p;
		};
    private:
        Params p;
    };
};
