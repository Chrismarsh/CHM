#include "water_flux.hpp"
#include "PhysConst.h"
namespace Units = PhysConst::units;
template<typename FluxType>
water_flux<FluxType>::water_flux(const double mm_per_s) : _value_mm_per_s(mm_per_s) {};

template<typename FluxType> water_flux<FluxType> water_flux<FluxType>::from_W_per_m_squared(const double Q,const Units::Celsius T)
{
    auto val_mm_per_s = Q / (PhysConst::water_reference_density() * PhysConst::Lv(T)) * mm_per_m;

    return water_flux(val_mm_per_s);    
};

template<typename FluxType> water_flux<FluxType> water_flux<FluxType>::from_mm_per_dt(const double WE,const size_t dt)
{
    // mm / step (seconds/step) -> mm / si
    auto val_mm_per_s = WE / dt;
    return water_flux(val_mm_per_s);
};

template<typename FluxType> water_flux<FluxType> water_flux<FluxType>::from_mm_per_s(const double WE)
{
    return from_mm_per_dt(WE,1);
};

template<typename FluxType> water_flux<FluxType> water_flux<FluxType>::from_m_per_s(const double WE)
{
    constexpr auto MM_PER_M = 1000.0;
    auto val_mm_per_s = WE * MM_PER_M;
    return from_mm_per_s(val_mm_per_s);
};

template<typename FluxType>
const double water_flux<FluxType>::mm_per_dt(const size_t dt) const
{
    return _value_mm_per_s * dt;
};

template<typename FluxType>
const double water_flux<FluxType>::mm_per_s() const
{
    return mm_per_dt(1.0);
};

template<typename FluxType>
const double water_flux<FluxType>::m_per_s() const
{
    return _value_mm_per_s / 1000.0;
};

template<typename FluxType>
const double water_flux<FluxType>::W_per_m_squared(const Celsius T) const
{
    PhysConst::units::Celsius temp{T};
    return _value_mm_per_s * PhysConst::water_reference_density() * PhysConst::Lv(temp) / mm_per_m;     
};

template<typename FluxType>
bool water_flux<FluxType>::empty() const { return _value_mm_per_s == 0.0; };
