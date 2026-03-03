#include <gtest/gtest.h>
#include "PhysConst.h"
#include "water_flux.hpp"

TEST(MassFluxTest,ConstructFromWperMSquared)
{
	double Q = 1400.0;
    water_flux m = water_flux<FluxType::latent>::from_W_per_m_squared(Q);
};

TEST(MassFluxTest,ReturnsMMperDT)
{
    double Q = 1400.0;
    PhysConst::units::Celsius T{0.0};
    water_flux m = water_flux<FluxType::latent>::from_W_per_m_squared(Q,T);

    auto result = m.W_per_m_squared(T);

    EXPECT_EQ(Q,result);
};

TEST(MassFluxTest,WperMSquaredtoMMperDT)
{
    auto Q = 250.5;
    PhysConst::units::Celsius T{10.0};
    water_flux m = water_flux<FluxType::latent>::from_W_per_m_squared(Q,T);

    int dt = 3600;
    double Q_massflux = m.mm_per_dt(dt);

    constexpr auto mm_per_m = 1000.0;

    // [Q] = J/m^2/s
    // [Q_massflux] = mm/s
    // [water density] = kg/m^3
    // [latent heat of vapourization] = J/kg
    // [milimetre per metre] = mm/m
    // [dt] = s/step
    //
    // J/m^2/s * m^3/kg * kg/J * mm/m * s/step =
    // J*m/s * 1/kg * kg/J * mm/m * s/step = 
    // m/s * mm/m * s/step = 
    // mm/step 
    auto Q_expect = Q / (PhysConst::water_reference_density() * PhysConst::Lv(T)) * mm_per_m;
    
    Q_expect *= dt;
    EXPECT_DOUBLE_EQ(Q_expect,Q_massflux);
};

TEST(MassFluxTest,GetBackWperMSquared)
{
    double Q = 1400.0;
    water_flux m = water_flux<FluxType::latent>::from_W_per_m_squared(Q);

    auto Q_expect = m.W_per_m_squared();

    EXPECT_DOUBLE_EQ(Q,Q_expect);
};

TEST(MassFluxTest,ConstructFromMMperS)
{
    double val = 200.0;
    water_flux m = water_flux<>::from_mm_per_s(val);
};

TEST(MassFluxTest,GetBackMMperS)
{
    const double val = 200.0;
    water_flux m = water_flux<FluxType::latent>::from_mm_per_s(val);

    auto result = m.mm_per_dt(1);                                  

    EXPECT_DOUBLE_EQ(result,val);

    auto result2 = m.mm_per_s(); // alias for mm_per_dt

    //result = result2 is implicit
    EXPECT_DOUBLE_EQ(result2,val);
};

TEST(MassFluxTest,ConstructFromMMperDT)
{
    auto val = 34.0;
    auto dt = 1000;
    water_flux m = water_flux<>::from_mm_per_dt(34.0,dt);

    auto result = m.mm_per_dt(dt);

    EXPECT_DOUBLE_EQ(result,val);

    auto result_per_s = m.mm_per_s();
    
    EXPECT_DOUBLE_EQ(val,result_per_s * dt);

};

TEST(MassFluxTest,ConstructFromMperS)
{
    auto input = 550.0;

    water_flux m = water_flux<>::from_m_per_s(input);

    auto result = m.mm_per_s();

    constexpr auto MM_PER_M = 1000.0;
    EXPECT_DOUBLE_EQ(result,input*MM_PER_M);
};

TEST(MassFluxTest,GetBackMperS)
{
    auto input = 550.0;

    water_flux m = water_flux<>::from_mm_per_s(input);

    auto result = m.m_per_s();

    constexpr auto M_PER_MM = 1 / 1000.0;
    EXPECT_DOUBLE_EQ(result,input * M_PER_MM);
};

TEST(MassFluxTest,PlusOperatorOverload)
{
	auto value1 = 50.4;
	auto m1 = water_flux<>::from_mm_per_s(value1);
	auto value2 = 250.3;
	auto m2 = water_flux<>::from_mm_per_s(value2);

	auto m3 = m1+m2;
	auto m4 = m2+m1;
	EXPECT_EQ(m3.mm_per_s(),value1+value2);
	EXPECT_EQ(m4.mm_per_s(),value1+value2);
};

TEST(MassFluxTest,MinusOperatorOverload)
{
	auto value1 = 50.4;
	auto m1 = water_flux<>::from_mm_per_s(value1);
	auto value2 = 250.3;
	auto m2 = water_flux<>::from_mm_per_s(value2);
	
	auto m3 = m1-m2;
	auto m4 = m2-m1;

	EXPECT_EQ(m3.mm_per_s(),value1 - value2);
	EXPECT_EQ(m4.mm_per_s(), value2 - value1);
};
