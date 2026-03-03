#include "katabatic_routing_glacier.hpp"
#include "PhysConst.h"

// data and view constructors
// module constructor
// module init
// module run
Units::Kelvin katabatic_routing_glacier::data::glacier_temperature()
{
	static double T = cfg_.get<double>("glacier_temp"_s,273.15);

	return Units::Kelvin{T};
};
Units::Kelvin katabatic_routing_glacier::data::air_temperature()
{
	update_value( [this]() -> auto& { return cache_->air_temperature; },
			[this]() { return (*face)["t"_s]; } );

	return Units::Kelvin{cache_->air_temperature};
};
Units::Pa katabatic_routing_glacier::data::air_pressure()
{
	update_value( [this]() -> auto& { return cache_->air_pressure; },
			[this]() { return (*face)["air_pressure"_s]; } );

	return Units::Pa{cache_->air_pressure};
};
Units::Pa katabatic_routing_glacier::data::vapour_pressure()
{
	update_value( [this]() -> auto& { return cache_->vapour_pressure; },
			[this]() { return (*face)["vapour_pressure"_s]; } );

	return Units::Pa{cache_->vapour_pressure};
};
Units::Pa katabatic_routing_glacier::data::vapour_pressure_surface()
{
	update_value( [this]() -> auto& { return cache_->vapour_pressure_surface; },
			[this]() { return (*face)["vapour_pressure_surface"_s]; } );

	return Units::Pa{cache_->vapour_pressure_surface};
};
Units::LapseRateSI katabatic_routing_glacier::data::lapse_rate()
{
	update_value( [this]() -> auto& { return cache_->lapse_rate; },
			[this]() { return (*face)["lapse_rate"_s]; } );

	return Units::LapseRateSI{cache_->lapse_rate};
};
void katabatic_routing_glacier::data::latent_heat(const double v)
{
	set_output( [this]() -> auto& { return cache_->latent_heat; },
			v);
};
void katabatic_routing_glacier::data::sensible_heat(const double v)
{
	set_output( [this]() -> auto& { return cache_->sensible_heat; },
			v);
};

const Units::Milimeters katabatic_routing_glacier::data::snowmelt()
{
	update_value( [this]() -> auto& { return cache_->snowmelt; },
			[this]() { return (*face)["snowmelt"_s]; } );

	return Units::Milimeters{cache_->snowmelt};
};
const Units::Milimeters katabatic_routing_glacier::data::firnmelt()
{
	if (!cache_ || std::isnan(cache_->firnmelt))
		return Units::Milimeters{0.0};
	
	return Units::Milimeters{cache_->firnmelt};
};
const Units::Milimeters katabatic_routing_glacier::data::icemelt()
{
	if (!cache_ || std::isnan(cache_->icemelt))
		return Units::Milimeters{0.0};
	
	return Units::Milimeters{cache_->icemelt};
};
void katabatic_routing_glacier::data::snowmelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->snowmelt_delayed; },
			v);
};
void katabatic_routing_glacier::data::firnmelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->firnmelt_delayed; },
			v);

};
void katabatic_routing_glacier::data::icemelt_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->icemelt_delayed; },
			v);
};
void katabatic_routing_glacier::data::total_delayed(double v)
{
	set_output( [this]() -> auto& { return cache_->total_delayed; },
			v);
};

const Units::Milimeters katabatic_routing_glacier::data::swe()
{
	update_value( [this]() -> auto& { return cache_->swe; },
			[this]() { return (*face)["swe"_s]; } );

	return Units::Milimeters{cache_->swe};
};
const Units::Watts_per_m2 katabatic_routing_glacier::data::melt_energy()
{
	auto Qsun = (*face)["Qsun"_s];
	auto Qrain = (*face)["Qrain"_s];
	
	if (!cache_)
		return Units::Watts_per_m2{Qsun + Qrain};

	auto Q_sensible = 0.0;	
	auto Q_latent = 0.0;

	if ( !std::isnan(cache_->latent_heat) )
		Q_latent = cache_->latent_heat;

	if ( !std::isnan(cache_->sensible_heat) )
		Q_sensible = cache_->sensible_heat;

	return Units::Watts_per_m2{Qsun + Qrain + Q_sensible + Q_latent};
};
bool katabatic_routing_glacier::data::update_now()
{
	size_t day = cfg_.get<size_t>("change_day"_s,300u);
	return day == global_param->day(); 
};

void katabatic_routing_glacier::data::glacier_water_equivalent(const double out)
{
	set_output( [this]() -> auto& { return cache_->glacier_water_equivalent; },
			out);
};
void katabatic_routing_glacier::data::total_depth(const double out)
{
	set_output( [this]() -> auto& { return cache_->total_depth; },
			out);
};
void katabatic_routing_glacier::data::firn_melt(const double out)
{
	set_output( [this]() -> auto& { return cache_->firnmelt; },
			out);
};
void katabatic_routing_glacier::data::ice_melt(const double out)
{
	set_output( [this]() -> auto& { return cache_->icemelt; },
			out);
};

