#include <gtest/gtest.h>
#include "Glacier.hpp"
#include "Atmosphere.h"

using namespace Glacier;

TEST(GlacierParams,ParamsAreConstExpr)
{
    static_assert(Params::g > 0.0);
    static_assert(Params::K > 0.0);
    static_assert(Atmosphere::psychrometric_constant > 0.0);
    static_assert(Atmosphere::water_density > 0.0);
    static_assert(Atmosphere::heat_capacity_air > 0.0);
    static_assert(Atmosphere::latent_heat_vapourization > 0.0);
    static_assert(Atmosphere::gas_constant_dry() > 0.0);
    static_assert(Atmosphere::gas_constant_vapour() > 0.0);
};

TEST_F(KatabaticParam,ParamsAreRunTime)
{
    auto& p = params();
    
    p.Ts_glacier = 300.0;
    p.Prandtl_number = 10.0;
     
    static_assert(Atmosphere::gas_constant_dry > 0.0);
    static_assert(Atmosphere::gas_constant_vapour > 0.0);
};

TEST_F(KatabaticParam,ParamsGettersAreMutableOrConst)
{
    /*
     *  const Params& and Params& getters to allow mutability to be explicit
     */ 
    static_assert(std::is_same_v<decltype(get_params()), const Params&>); 
    static_assert(std::is_same_v<decltype(get_mutable_params()), Params&>); 

    // Test that changing the mutable version affects the immutable version
    constexpr auto T = 200;
    {
        auto& mut_p = get_mutable_params();

        mut_p.glacier_surface_temp = T;
    }
    {
        auto& p = get_params();
    }
    
    EXPECT_EQ(p.glacier_surface_temp,T);
};

class KatabaticParam : public ::testing::Test
{
protected:
    katabatic_melt_energy katabatic;
    inline constexpr auto air_temperature = -5.0;
    inline constexpr auto lapse_rate = 10.0;
};

TEST_F(KatabaticParam,ReturnsNotZero)
{ 
    State::defaults();
    auto K = katabatic.K(air_temperature,lapse_rate);

    EXPECT_NE(K,0.0);
}

TEST_F(KatabaticParam,AirTempEqualGlacierTempReturnsZero)
{
    auto& p = get_params();
    auto K = katabatic.K(p.glacier_surface_temp,lapse_rate);

    EXPECT_EQ(K,0.0);
};

TEST_F(KatabaticParam,NegativeLapseReturnsNaN)
{
    auto K = katabatic.K(air_temperature,-lapse_rate);

    EXPECT_TRUE(std::isnan(K));
};

class KatabaticLatent : public KatabaticParam
{
protected:
    inline constexpr auto air_pressure = 101.325;
    inline constexpr auto vapour_pressure = air_pressure * 0.4;
};

TEST_F(KatabaticLatent,ReturnsNonZero)
{
    auto Q = katabatic.Q_latent(air_pressure,vapour_pressure,air_temperature);

    EXPECT_NE(Q,0.0); 
};

TEST_F(KatabaticLatent,AirTempEqualGlacierTempReturnsZero)
{
    auto& p = get_params();

    auto Q = katabatic.Q_latent(air_pressure,vapour_pressure,p.glacier_surface_temp);

    EXPECT_EQ(Q,0.0);
};

TEST_F(KatabaticLatent,PequalVapourPressureReturnsZero)
{
    auto Q = katabatic.Q_latent(air_pressure,air_pressure,air_temperature);

    EXPECT_EQ(Q,0.0);
};

class KatabaticSensible : public KatabaticLatent
{
};

TEST_F(KatabaticSensible,ReturnsNonZero)
{
    auto Q = katabatic.Q_sensible(vapour_pressure);

    EXPECT_NE(Q,0.0);
};

TEST_F(KatabaticSensible,SurfaceEqualsAirWaterVapourVapourReturnsZero)
{
    auto& p  = get_params();

    auto Q = katabatic.Q_sensible(p.water_vapour_pressure);

    EXPECT_EQ(Q,0.0);
};

class KatabaticRainEnergy : public KatabaticParam
{

};

TEST_F(KatabaticRainEnergy,ReturnsNonZero)
{
    auto Q = katabatic.Q_precip(rainfall,rain_temperature);

    EXPECT_NE(Q,0.0);
};

TEST_F(KatabaticRainEnergy,RainTempEqualGlacierTempReturnsZero)
{
    auto& p = get_params();
    auto Q = Katabatic.Q_precip(rainfall,p.glacier_surface_temp);

    EXPECT_EQ(Q,0.0);
};

TEST_F(KatabaticRainEnergy,RainTempEqualFreezingTempReturnsZero)
{
    auto Q = Katabatic.Q_precip(rainfall,0.0);

    EXPECT_EQ(Q,0.0);
};
