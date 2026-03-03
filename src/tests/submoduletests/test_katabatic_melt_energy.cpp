#include <gtest/gtest.h>
#include "katabatic_melt_energy.hpp"
#include "PhysConst.h"
#include "Atmosphere.h"
#include <array>
#include <stdexcept>

using namespace katabatic_melt_energy;
namespace Units = PhysConst::units;

inline Params test_params()
{
	Params p;
	
	p.prandtl = 5.0;
	p.k = 4e-4;
	p.k2 = 1.0;
	p.seconds_per_step = 3600;

	return p;
};

TEST(ParamsKatabaticTest,ParamsSetByPhysConst)
{
    Params p;


    EXPECT_EQ(p.heat_capacity_air,PhysConst::Cp());
    EXPECT_EQ(p.molecular_wt_ratio,PhysConst::M());
    EXPECT_EQ(p.g,PhysConst::g());
}

class AirDensityTest : public ::testing::Test
{
protected:
    struct Data
    {
        Units::Pa pressure;
        Units::Pa vapour_pressure;
        Units::Kelvin air_temperature;
        Units::DensitySI expected_air_density;

        Data(const double p, const double e, const double T, const double rho) : 
            pressure{p}, vapour_pressure{e}, air_temperature{T}, expected_air_density{rho} {};
    };

    std::array<Data,2> data = {{
        {1e5,1e4,273.15,1.227141074763899},
        {1.1e5,1e3,300.2,1.2720886528749895}
    }};
};
TEST_F(AirDensityTest,Analytical)
{
    for (auto d : data)
    {
        auto result = Atmosphere::air_density(d.pressure,d.air_temperature,d.vapour_pressure);
        EXPECT_DOUBLE_EQ(result.value,d.expected_air_density.value);
    };

};

TEST_F(AirDensityTest,ZeroPressureThrows)
{
    EXPECT_THROW(Atmosphere::air_density( Units::Pa{0.0},
                Units::Kelvin{273.15},
                Units::Pa{1e7}),std::runtime_error);
};

class BulkCoefficientTest : public ::testing::Test
{
protected:
    Params p = test_params();

    struct Input
    {
        double C;
        double gamma;
        double prandtl;
        double glacier_temperature;
        double solution;
    };
    std::array<Input,2> inputs{{
        {-10,0.005,5.0,273.15,0.004791841256996453},
        {10.0,0.001,1.0,250.3,-0.025028948481418418}
    }};
};

TEST_F(BulkCoefficientTest,AnalyticalBulkCoefficent)
{
    for (auto input : inputs)
    {
		input.solution = p.k * std::pow(p.k2,2.0) * input.C * std::sqrt(p.g / (input.glacier_temperature * input.gamma * input.prandtl));
		p.prandtl = input.prandtl;
        auto result = bulk_coefficient( p,
                Units::TempDiff{input.C},
                Units::LapseRateSI{input.gamma},
                Units::Kelvin{input.glacier_temperature});

        EXPECT_EQ(result.m_per_s(),input.solution);
    }

};

TEST_F(BulkCoefficientTest,CertainZerosThrow)
{
    /*
     * inputs[0] -> gamma  = 0, divide by zero
     * inputs[1] -> Prandtl number = 0, divide by zero
     * inputs[2] -> glacier temperature = 0, divide by zero
     */
    std::array<Input,3> inputs {{
        {1.0,0.0,1.0,200.0,0.0}, 
        {1.0,1.0,0.0,200.0,0.0}, 
        {1.0,1.0,1.0,0.0,0.0}
    }};

    for (auto input : inputs)
    {
		p.prandtl = 0.0;
        EXPECT_THROW(bulk_coefficient( p,
                Units::Kelvin{input.C},
                Units::LapseRateSI{input.gamma},
                Units::Kelvin{input.glacier_temperature});,
                std::runtime_error);
    }

};

class KatabaticHeatTest : public ::testing::Test
{
protected:
    Params p = test_params(); 

};

TEST_F(KatabaticHeatTest,AnalyticalSensibleHeat)
{
    const auto coefficient = water_flux<>::from_m_per_s(1e3);
    const Units::Celsius deficit{10.0};
    const Units::DensitySI air_density{1.5};
    const auto expected = air_density.value * p.heat_capacity_air * coefficient.m_per_s() * deficit.value;

    auto result = sensible_heat(p,coefficient,deficit,air_density);

    EXPECT_DOUBLE_EQ(result.W_per_m_squared(),expected);
    
};

TEST_F(KatabaticHeatTest,AnalyticlLatentHeat)
{
    const auto coefficient = water_flux<>::from_m_per_s(45.0);
    const Units::DensitySI air_density{1.102};
    const Units::Celsius glacier_temperature{-1.3};
    const Units::Pa vapour_pressure_deficit{1e-3};
    const auto expected = p.molecular_wt_ratio * air_density.value * PhysConst::Lv(glacier_temperature) * coefficient.m_per_s() * vapour_pressure_deficit.value;

    auto result = latent_heat(p,coefficient,vapour_pressure_deficit,glacier_temperature,air_density);

    EXPECT_EQ(result.W_per_m_squared(),expected);
};


// Simple mock data class that satisfies KatabaticData concept
class MockKatabaticData {
public:
    // Store the values that will be passed to latent_heat and sensible_heat
    double last_latent_heat_value{0};
    double last_sensible_heat_value{0};
    
    // Required getters
    Units::Kelvin glacier_temperature() const { return glacier_temp; }
    Units::Kelvin air_temperature() const { return air_temp; }
    Units::Pa air_pressure() const { return air_press; }
    Units::Pa vapour_pressure() const { return vapour_press; }
    Units::Pa vapour_pressure_surface() const { return vapour_press_surface; }
    Units::LapseRateSI lapse_rate() const { return lapse_rate_val; }
    
    // Methods called by Model
    void latent_heat(double value) { last_latent_heat_value = value; }
    void sensible_heat(double value) { last_sensible_heat_value = value; }
    
    // Setters for test configuration
    void set_glacier_temperature(Units::Kelvin v) { glacier_temp = v; }
    void set_air_temperature(Units::Kelvin v) { air_temp = v; }
    void set_air_pressure(Units::Pa v) { air_press = v; }
    void set_vapour_pressure(Units::Pa v) { vapour_press = v; }
    void set_vapour_pressure_surface(Units::Pa v) { vapour_press_surface = v; }
    void set_lapse_rate(Units::LapseRateSI v) { lapse_rate_val = v; }
    
private:
    Units::Kelvin glacier_temp{273.15};
    Units::Kelvin air_temp{273.15};
    Units::Pa air_press{101325};
    Units::Pa vapour_press{1000};
    Units::Pa vapour_press_surface{500};
    Units::LapseRateSI lapse_rate_val{0.0065};
};

class KatabaticModelTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up default parameters
        params.seconds_per_step = 3600;  // 1 hour
    }
    
    Params params;
};

TEST_F(KatabaticModelTest, ExecuteImplComputesCorrectSensibleHeat) {
    MockKatabaticData data;
    data.set_air_temperature(Units::Kelvin{275.15});  // 2°C
    data.set_glacier_temperature(Units::Kelvin{273.15});  // 0°C
    data.set_air_pressure(Units::Pa{101325});
    data.set_vapour_pressure(Units::Pa{1000});
    data.set_lapse_rate(Units::LapseRateSI{0.0065});
    
    Model<MockKatabaticData> model;
    
    // Manually compute expected values using the static functions
    auto rho_air = Atmosphere::air_density(
        data.air_pressure(),
        data.glacier_temperature(),
        data.vapour_pressure()
    );
    
    Units::Kelvin T_deficit{
        data.air_temperature().value - data.glacier_temperature().value
    };
    
    auto K = bulk_coefficient(
        params,
        T_deficit,
        data.lapse_rate(),
        data.glacier_temperature()
    );
    
    auto expected_sensible = sensible_heat(params, K, T_deficit, rho_air);
    auto expected_sensible_mm = expected_sensible.mm_per_dt(params.seconds_per_step);
    
    // Execute the model
    model.execute_impl(data);
    
    // Verify
    EXPECT_DOUBLE_EQ(data.last_sensible_heat_value, expected_sensible_mm);
}

TEST_F(KatabaticModelTest, ExecuteImplComputesCorrectLatentHeat) {
    MockKatabaticData data;
    data.set_air_temperature(Units::Kelvin{275.15});
    data.set_glacier_temperature(Units::Kelvin{273.15});
    data.set_air_pressure(Units::Pa{101325});
    data.set_vapour_pressure(Units::Pa{1000});
    data.set_vapour_pressure_surface(Units::Pa{500});
    data.set_lapse_rate(Units::LapseRateSI{0.0065});
    
    Model<MockKatabaticData> model;
    
    // Manually compute expected values
    auto rho_air = Atmosphere::air_density(
        data.air_pressure(),
        data.glacier_temperature(),
        data.vapour_pressure()
    );
    
    Units::TempDiff T_deficit{
        data.air_temperature().value - data.glacier_temperature().value
    };
    
    auto K = bulk_coefficient(
        params,
        T_deficit,
        data.lapse_rate(),
        data.glacier_temperature()
    );
    
    Units::Pa vapour_pressure_deficit{
        data.vapour_pressure().value - data.vapour_pressure_surface().value
    };
    
    auto expected_latent = latent_heat(
        params,
        K,
        vapour_pressure_deficit,
        data.glacier_temperature(),
        rho_air
    );
    auto expected_latent_mm = expected_latent.mm_per_dt(params.seconds_per_step);
    
    // Execute
    model.execute_impl(data);
    
    // Verify
    EXPECT_DOUBLE_EQ(data.last_latent_heat_value, expected_latent_mm);
}

TEST_F(KatabaticModelTest, ExecuteImplWithZeroTemperatureGradient) {
    MockKatabaticData data;
    data.set_air_temperature(Units::Kelvin{273.15});  // Same as glacier
    data.set_glacier_temperature(Units::Kelvin{273.15});
    data.set_air_pressure(Units::Pa{101325});
    data.set_vapour_pressure(Units::Pa{1000});
    
    Model<MockKatabaticData> model;
    model.execute_impl(data);
    
    // With zero temperature gradient, sensible heat should be zero
    EXPECT_DOUBLE_EQ(data.last_sensible_heat_value, 0.0);
}

TEST_F(KatabaticModelTest, ExecuteImplWithZeroVapourPressureDeficit) {
    MockKatabaticData data;
    data.set_vapour_pressure(Units::Pa{1000});
    data.set_vapour_pressure_surface(Units::Pa{1000});  // Equal pressures
    data.set_air_temperature(Units::Kelvin{275.15});
    data.set_glacier_temperature(Units::Kelvin{273.15});
    
    Model<MockKatabaticData> model;
    model.execute_impl(data);
    
    // With zero vapour pressure deficit, latent heat should be zero
    EXPECT_DOUBLE_EQ(data.last_latent_heat_value, 0.0);
}
