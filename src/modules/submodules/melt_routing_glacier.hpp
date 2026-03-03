#pragma once

#include "base_step.hpp"
#include "PhysConst.h"
#include <concepts>

namespace Units = PhysConst::units;

namespace GlacierRouting
{
	struct Params {
		size_t chain_length = 10u;
		double k = 4e3;
		size_t seconds_per_step = 3600u;
	};

	class LinearReservoir {
		static const double get_c0(const double k, const size_t dt);
		static const double get_c1(const double k, const size_t dt);
		const double c0;
		const double c1;

		double last_input;
		double last_output;
	public:
		LinearReservoir(const double k, const double dt);

		const double step(const double input);
	};

	class GlacierReservoir
	{
		const Params* const p;
		std::vector<LinearReservoir> reservoir_chain;
		std::vector<double> old_outputs;
	public:
		GlacierReservoir(const Params* _p);

		const double step(const double input);
	};

	struct State 
	{
		State(const Params* p);
		// Delete copy constructors
		State(const State&) = delete;
		State& operator=(const State&) = delete;
		GlacierReservoir snow_route;
		GlacierReservoir firn_route;
		GlacierReservoir ice_route;
	};

	template<class T>
	concept GlacierRoutingData = requires(T& t)
	{
		{ t.get_state() } -> std::same_as<State&>;
		{ t.snowmelt() } -> std::same_as<const Units::Milimeters>;
		{ t.firnmelt() } -> std::same_as<const Units::Milimeters>;
		{ t.icemelt() } -> std::same_as<const Units::Milimeters>;

		{ t.snowmelt_delayed(std::declval<double>()) } -> std::same_as<void>;
		{ t.firnmelt_delayed(std::declval<double>()) } -> std::same_as<void>;
		{ t.icemelt_delayed(std::declval<double>()) } -> std::same_as<void>;
		{ t.total_delayed(std::declval<double>()) } -> std::same_as<void>;
	};

	template<GlacierRoutingData Data>
	class Model : public base_step<Model<Data>,Data>
	{
		Params p;
	public:
		void execute_impl(Data& d) const
		{
			auto snowmelt = d.snowmelt();
			auto firnmelt = d.firnmelt();
			auto icemelt = d.icemelt();
			auto& s = d.get_state();
			
			auto snow = s.snow_route.step(snowmelt.value);
			auto firn = s.firn_route.step(firnmelt.value);
			auto ice = s.ice_route.step(icemelt.value);

			d.snowmelt_delayed(snow);
			d.firnmelt_delayed(firn);
			d.icemelt_delayed(ice);

			d.total_delayed(snow + firn + ice);
		};

		Params& get_params() { return p; }
	};
}
