#pragma once
#include "PhysConst.h"
#include <cstddef>

namespace Units = PhysConst::units;
enum class FluxType {
latent,
sensible,
Default
};

template<FluxType type = FluxType::Default>
class water_flux
{
	static_assert(!std::constructible_from<water_flux> &&
            !std::constructible_from<water_flux,double>);

    double _value_mm_per_s;
    explicit water_flux(const double mm_per_s);

    static constexpr auto mm_per_m = 1000.0;
    
    using Celsius = PhysConst::units::Celsius;
public:

    // Defaults to freezing temperature
    static water_flux from_W_per_m_squared(const double Q, 
            const Celsius T = Celsius{0.0} );
    static water_flux from_mm_per_dt(const double WE,const size_t dt);    
    static water_flux from_mm_per_s(const double WE);
    static water_flux from_m_per_s(const double WE);

    const double mm_per_dt(const size_t dt) const ;
    const double W_per_m_squared(const Celsius T = Celsius{0.0}) const 
		requires (type == FluxType::latent);//std::is_same_v<FluxType,flux::latent>;

    const double mm_per_s() const;
    const double m_per_s() const;
    bool empty() const;

	water_flux operator+(const water_flux& other) const
	{
		return water_flux(_value_mm_per_s + other._value_mm_per_s);
	};

	water_flux operator-(const water_flux& other) const
	{
		return water_flux(_value_mm_per_s - other._value_mm_per_s);
	};

	void operator+=(const water_flux& other)
	{
		_value_mm_per_s += other._value_mm_per_s;
	};

	void operator-=(const water_flux& other)
	{
		_value_mm_per_s -= other._value_mm_per_s;
	};
};

template<FluxType type>
water_flux<type>::water_flux(const double mm_per_s) : _value_mm_per_s(mm_per_s) {};

template<FluxType type> 
water_flux<type> water_flux<type>::from_W_per_m_squared(const double Q)
{
	auto val_mm_per_s = 0.0;
	if constexpr (type == FluxType::latent)
		val_mm_per_s = Q / (PhysConst::water_reference_density() * PhysConst::Lf()) * mm_per_m;
	else 
	{
		static_assert(type != FluxType::latent && 
				"W/m^2 functions are not supported for sensible heat or other quantities"
				"If a sensible heat conversion to mm/dt is required, it must be done manually"
				"Carefully consider if this path is physical as sensible heat == temperature change, not a mass flux.");
	}

    return water_flux(val_mm_per_s);    
};

template<FluxType type> 
water_flux<type> water_flux<type>::from_mm_per_dt(const double WE,const size_t dt)
{
    // mm / step (seconds/step) -> mm / si
    auto val_mm_per_s = WE / dt;
    return water_flux(val_mm_per_s);
};

template<FluxType type> 
water_flux<type> water_flux<type>::from_mm_per_s(const double WE)
{
    return from_mm_per_dt(WE,1);
};

template<FluxType type> 
water_flux<type> water_flux<type>::from_m_per_s(const double WE)
{
    constexpr auto MM_PER_M = 1000.0;
    auto val_mm_per_s = WE * MM_PER_M;
    return from_mm_per_s(val_mm_per_s);
};

template<FluxType type>
const double water_flux<type>::mm_per_dt(const size_t dt) const
{
    return _value_mm_per_s * dt;
};

template<FluxType type>
const double water_flux<type>::mm_per_s() const
{
    return mm_per_dt(1.0);
};

template<FluxType type>
const double water_flux<type>::m_per_s() const
{
    return _value_mm_per_s / 1000.0;
};

template<FluxType type>
const double water_flux<type>::W_per_m_squared(const Celsius T) const
requires (type == FluxType::latent)
{
    PhysConst::units::Celsius temp{T};
    return _value_mm_per_s * PhysConst::water_reference_density() * PhysConst::Lv(temp) / mm_per_m;     
};

template<FluxType type>
bool water_flux<type>::empty() const { return _value_mm_per_s == 0.0; };
