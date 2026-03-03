#include <concepts>
#include <gtest/gtest.h>
#include <memory>
#include "Glacier.hpp"
#include "PhysConst.h"
#include <memory>
#include <numeric>
#include <stdexcept>

using namespace Glacier;
inline constexpr auto MARGIN = 1e-12;
/*
 * ParamsTest: Verifies the Params struct initialization.
 */

TEST(ParamsTest,ParamsSetByPhysConst)
{
    Params p;

	EXPECT_EQ(p.water_density,PhysConst::water_reference_density());
	EXPECT_EQ(p.heat_capacity_air,PhysConst::Cp());
	EXPECT_EQ(p.latent_heat_vapour,PhysConst::Lv());
	EXPECT_EQ(p.gas_constant_dry,PhysConst::RgasDry());
	EXPECT_EQ(p.gas_constant_vapour,PhysConst::RgasVapour());
	EXPECT_EQ(p.molecular_wt_ratio,PhysConst::M());
    
    auto _ = p.densify_version;
};

/*
 * This function exists so that changes to the default values in the 
 * Glacier::Params struct can be changed without impacting the test success.
 */ 
[[maybe_unused]] static inline Params test_default_params()
{
    Params p;
    p.densify_version  = DensifyVersion::HerronLangway;
    p.small_increment = 25.0;
    p.big_increment = 50.0;
    p.critical_density = 550.0;
	p.firn_to_ice_density = 830.0;
	p.seconds_per_step = 3600.0;
    return p;
};

/*
 *  DensificationTest: Verifies the calculations of densification of firn.
 *
 *  Two methods: Linear method, and method based on Herron and Langway (1980)
 */ 

class DensificationTest : public ::testing::Test
{
protected:
    Params p;

    const double critical_density = 550.0;
    
};

/*
 * DensificationTest - HerronLangway: Tests the Herron and Langway (1980) firn
 * densification algorithm.
 * 
 */
class HerronLangwayTest : public DensificationTest
{
    static_assert(PhysConst::rho_ice() == 917.0, 
            "The original model assumes this specific value for ice density and cannot be modified.\n"
            "Consider refactoring if PhysConst::rho_ice() must be modified.\n"
            "I would suggest creating a ice_density constexpr function in Glacier.hpp in the Glacier namespace.");
};

TEST_F(HerronLangwayTest,ValidationAgainstPrecalculated)
{
    std::array<double,4> heights{1.0,1.2,0.9,0.5};
    std::array<double,4> densities{250.0,320.0,400.0,810.0};
    std::array<double,4> glacier_temperature{273.15,250.7,250.8,200.1};
    std::array<double,4> Accumulate_rate{0.11,0.265,0.4,0.022};
    const std::array<double,4> new_density{
        271.4597957366285,
        356.16647094729564,
        454.56304130714364,
        688.3589090630271
    };
    
    double depth = 0.0;
    for (size_t layer = 0; layer < densities.size(); ++layer)
    {
        depth += heights[layer];
        auto result = Densification::HerronLangway(depth,densities[layer],glacier_temperature[layer],Accumulate_rate[layer],critical_density);

        EXPECT_DOUBLE_EQ(result,new_density[layer]);
    };

};

TEST_F(HerronLangwayTest,MaxDensityThrows)
{
    double height = 1.0;
    double density = PhysConst::rho_ice() * 1.1;

    // density >= rho_ice = 917.0 kg/m^3 should not be permitted.
    EXPECT_THROW(
            Densification::HerronLangway(height,density,273.15,1.0,critical_density),
            std::runtime_error);
    EXPECT_THROW(
            Densification::HerronLangway(height,PhysConst::rho_ice(),273.15,1.0,critical_density),
            std::runtime_error);

};

/*
 * DensificationTest - Linear: Validates the linear densification algorithm used
 * for testing and baseline comparisons.
 * 
 */

TEST_F(DensificationTest,Linear)
{
    std::array<double,4> densities{250.0,400.0,725.0};
    const double small_increment = 25.0;
    const double big_increment = 55.0;

    auto result = Densification::Linear(densities[0],small_increment);
    EXPECT_EQ(result,densities[0] + small_increment);
    
    result = Densification::Linear(densities[1],small_increment);
    EXPECT_EQ(result,densities[1] + small_increment);

    result = Densification::Linear(densities[2],big_increment);
    EXPECT_EQ(result,densities[2] + big_increment);
};

/*/
 * FirnLayerTest: Class containing the data nad behaviour of a single layer of firn 
 */

class FirnLayerTest : public ::testing::Test
{
protected:
    // Testing on a single layer, constructed in each test
    std::unique_ptr<Layer> layer; 
    Params p = test_default_params();

    // placeholder to be allocated later and passed to the layer object
    Units::Milimeters water_equivalent;

    // Height is metres, density is kg/m^3
    static_assert(std::derived_from<Layer::Height,Units::Metres>);
    static_assert(std::derived_from<Layer::Density,Units::DensitySI>);
};

// Layer constructs from Units::Milimeters struct and returns the same value with the same type
TEST_F(FirnLayerTest,ConstructsFromMilimetreStruct)
{
    water_equivalent.value = 30.0;
    layer = std::make_unique<Layer>(water_equivalent);

    auto WE = layer->water_equivalent();

    EXPECT_EQ(WE.value,water_equivalent.value);
    
};

// Construct from Layer::Height and Layer::Density objects.
// Reproduce WE using WE formula: Dingman (2008), equation (5-12)
TEST_F(FirnLayerTest,ConstructsFromHeightandDensity)
{
    Layer::Height height{1.0};
    Layer::Density density{300.0};
    layer = std::make_unique<Layer>(height,density);

    auto WE = layer->water_equivalent();

    // returns water equivalent from original height/density
    auto constexpr MM_PER_M = 1000.0;
    EXPECT_EQ(WE.value,height.value * density.value / p.water_density * MM_PER_M);
};

TEST_F(FirnLayerTest,PartialRemove)
{
    water_equivalent.value = 350.0;
    layer = std::make_unique<Layer>(water_equivalent);

    Units::Milimeters to_remove{150.0};
    layer->remove(to_remove);

	// Layer uses a bisection method upon construction
	// to convert the WE to height and density. There will
	// definitely be a loss of precision.
	//
	// adding a layer is typically only once a year, so its acceptable
	// uncertainty.
    EXPECT_NEAR(layer->water_equivalent().value,water_equivalent.value - to_remove.value,MARGIN);
};

// Remove function is just for reducing a layer by an amount and should not be used for
// removal, thats what pop_back() is for!
TEST_F(FirnLayerTest,ExcessiveRemoveThrowsException)
{
    water_equivalent.value = 250;
    Units::Milimeters to_remove{water_equivalent.value * 2.0};
    layer = std::make_unique<Layer>(water_equivalent);

    EXPECT_THROW(layer->remove(to_remove),std::invalid_argument);
};

// Same as above, shouldn't be equal to total either. 
// A layer with no WE should be deleted as an object, not just set to zero
TEST_F(FirnLayerTest,RemoveEqualToTotalThrows)
{
    water_equivalent.value = 903.0;
    Units::Milimeters to_remove{water_equivalent.value};
    to_remove += water_equivalent;
    layer = std::make_unique<Layer>(water_equivalent);
    
    EXPECT_THROW(layer->remove(to_remove),std::invalid_argument);
};

// 
TEST_F(FirnLayerTest,LinearDensification)
{
    water_equivalent.value = 50.0;
    layer = std::make_unique<Layer>(water_equivalent);

    layer->densifyLinear(Layer::Density{p.small_increment},Layer::Density{p.big_increment},Layer::Density{p.critical_density});

    EXPECT_EQ(layer->water_equivalent().value,water_equivalent.value);

};

TEST_F(FirnLayerTest,HerronLangwayDensification)
{
    water_equivalent.value = 500.0;
    layer = std::make_unique<Layer>(water_equivalent);

    Units::Kelvin T{273.15};
    layer->densifyHerronLangway(T,0.0,1.2,Layer::Density{550.0});

    EXPECT_DOUBLE_EQ(layer->water_equivalent().value,water_equivalent.value);
};

/*/
 * LayeredFirnTest: Test fixture for the LayeredFirn class, which manages a stratified column
 * of firn layers.
 * 
 * In the language of C++ containers: the array "front" is the bottom of the firn and the back is the top of the firn.
 *
 * Uses a double-ended queue (deque) so that removal and adding to front/back are equally efficient.
 * deque is more efficient at adding new elements to the ends than std::vector. deque does not 
 * require that data is layed out sequentially, but vector does. Therefore, adding or removing a firn
 * layer from the front/back doesn't require moving the vector. It simply frees the bytes. 
 *
 * There are few loops over the layers, so it should be ok. 
 * 
 */
class LayeredFirnTest : public ::testing::Test
{
protected:
    Params p = test_default_params();
    std::unique_ptr<LayeredFirn> firn; 


    static constexpr auto MM_PER_M = 1000.0;

    std::vector<std::pair<Layer::Height,Layer::Density>> init_h_rho = {
            {Layer::Height{0.6},Layer::Density{250.0}},
            {Layer::Height{0.3},Layer::Density{350.0}},
            {Layer::Height{0.2},Layer::Density{840.0}},
            {Layer::Height{0.1},Layer::Density{900.0}},
            };


    std::deque<Layer> old_layers()
    {
        std::deque<Layer> result;
        
        auto expected_firn = 0.0; 
        
        for (auto init : init_h_rho)
        {
            result.push_front(Layer(init.first,init.second));
        };

        return result;
    };
};

TEST_F(LayeredFirnTest,DefaultConstructNoLayeredFirn)
{
    firn = std::make_unique<LayeredFirn>(&p);

    EXPECT_EQ(firn->water_equivalent().value,0.0);
};

TEST_F(LayeredFirnTest,ConstructFromLayersVector)
{
    std::vector<Layer> init_layers = [this]() { 
        auto layers = old_layers();
        std::vector<Layer> result(layers.begin(),layers.end());
        return result;
    }();

    auto expected_firn = 0.0; 
    
    for (auto init : init_h_rho)
    {
        expected_firn += init.first.value * init.second.value / PhysConst::water_reference_density() * MM_PER_M;
    }
    LayeredFirn firn2(&p,init_layers);
    firn = std::make_unique<LayeredFirn>(&p,init_layers);

    EXPECT_EQ(firn->water_equivalent().value,expected_firn);
    
};

TEST_F(LayeredFirnTest,Accumulate)
{
    constexpr Units::Milimeters WE{350.0}; 

    firn = std::make_unique<LayeredFirn>(&p);

    firn->accumulate(WE);

    auto result = firn->water_equivalent();

    EXPECT_DOUBLE_EQ(result.value,WE.value);

};

TEST_F(LayeredFirnTest,ConvertToIce)
{
    std::deque<Layer> layers = old_layers();

    firn = std::make_unique<LayeredFirn>(&p,layers); 
    Units::Milimeters init_WE = firn->water_equivalent();

    auto result = firn->convert_to_ice();

    auto& new_layers = firn->get_layers();

    EXPECT_EQ(firn->water_equivalent().value + result->value,
       init_WE.value); 

    EXPECT_EQ(layers.size() - 1,new_layers.size());

};

TEST_F(LayeredFirnTest,NoConvertToIce)
{
    std::deque<Layer> layers = old_layers();
    std::reverse(layers.begin(),layers.end());

    firn = std::make_unique<LayeredFirn>(&p,layers);
    Units::Milimeters init_WE = firn->water_equivalent();

    auto result = firn->convert_to_ice();

    auto& new_layers = firn->get_layers();

    EXPECT_EQ(firn->water_equivalent().value,init_WE.value);

    EXPECT_EQ(layers.size(),new_layers.size());

    EXPECT_EQ(layers.front().water_equivalent().value,
            new_layers.front().water_equivalent().value);

};


/*/
 * LayeredFirnMeltTest: Specialized test fixture for testing melt calculations of firn for the
 * variety of pathways 
 */

class LayeredFirnMeltTest : public LayeredFirnTest
{
protected:
	static constexpr std::array WE{200.0,150.0,303.0,525.0,201.0,25.0};	
	std::deque<Layer> layers;
	
	void SetUp() override {
		constexpr size_t N = WE.size();

		for (size_t n = 0; n < N; ++n)
		{
			Units::Milimeters water_eq{WE[n]};
			Layer layer(water_eq);
			layers.push_back(layer);
		};
		firn = std::make_unique<LayeredFirn>(&p,layers);
	};
};

/*
 * Unrelated to melt, but built deque here.
 */
TEST_F(LayeredFirnMeltTest,ConstructFromDeque)
{
    const auto sum = std::accumulate(WE.begin(),WE.end(),0.0);

    EXPECT_DOUBLE_EQ(sum,firn->water_equivalent().value);
    
    size_t ind = 0;
    for (const auto layer : firn->get_layers())
    {
        EXPECT_DOUBLE_EQ(layer.water_equivalent().value, WE[ind]);
        ++ind;
    }
};


TEST_F(LayeredFirnMeltTest,PartialMeltSingleLayer)
{
    firn = std::make_unique<LayeredFirn>(&p);
    constexpr Units::Milimeters input{150.0};
    constexpr double melt_flux{75.0};
    firn->accumulate(input);
    
    // second argument is dt in units of seconds / step
    // Must use p.seconds_per_step because internals do the same
    auto F = water_flux<FluxType::latent>::from_mm_per_dt(melt_flux,p.seconds_per_step);

    auto melt_info = firn->melt(F);

    EXPECT_DOUBLE_EQ(melt_info.melt.value,melt_flux);

    EXPECT_EQ(melt_info.remaining_energy.mm_per_dt(p.seconds_per_step),0.0);
    
    auto& layers = firn->get_layers();

    EXPECT_EQ(layers.size(),1);
    EXPECT_EQ(layers.front().water_equivalent(),layers.back().water_equivalent());
    EXPECT_NEAR(layers.at(0).water_equivalent().value,input.value - melt_info.melt.value,MARGIN);
};

TEST_F(LayeredFirnMeltTest,PartialMeltManyLayer)
{
	constexpr auto melt_energy_mm_per_dt = 705.0;
	water_flux m = water_flux<FluxType::latent>::
		from_mm_per_dt(melt_energy_mm_per_dt,p.seconds_per_step);
	double initial_WE = std::accumulate(WE.begin(),WE.end(),0.0);
	double expected_final_WE = initial_WE - melt_energy_mm_per_dt;
	
	MeltInfo info = firn->melt(m);

    auto sum  = 0.0;
    auto excess = 0.0;
    for (auto it = WE.rbegin(); it != WE.rend(); ++it)
    {
        sum += *it;
        if (sum >= melt_energy_mm_per_dt)
        {
            excess = sum - melt_energy_mm_per_dt;
            break;
        }
    }

    // amount melted is equal to input flux
	EXPECT_DOUBLE_EQ(info.melt.value,m.mm_per_dt(p.seconds_per_step));
	
    // final WE is correct
	EXPECT_DOUBLE_EQ(firn->water_equivalent().value,
			expected_final_WE);
    
    // layers removed at END only.
    // This isn't quite right. It will fail at the last iteration
    // num_layers()? Might also be useful to "predict" number of layers.
    const auto& layers = firn->get_layers();
    auto layers_it = layers.begin();
    auto WE_it = WE.begin();

    while(layers_it != layers.end()-1 && WE_it != WE.end())
    {
        EXPECT_DOUBLE_EQ((*layers_it).water_equivalent().value,*WE_it);
        ++layers_it;
        ++WE_it;
    }
    
    EXPECT_NEAR((*layers_it).water_equivalent().value,excess,MARGIN);
    ++layers_it;
    EXPECT_TRUE(layers_it == layers.end()); 
    EXPECT_FALSE(WE_it == WE.end());

};

TEST_F(LayeredFirnMeltTest,MeltAllLayersExcess)
{
	constexpr double initial_WE = std::accumulate(WE.begin(),WE.end(),0.0);
	constexpr auto melt_energy_mm_per_dt = initial_WE * 1.25;
	water_flux m = water_flux<FluxType::latent>::
		from_mm_per_dt(melt_energy_mm_per_dt,p.seconds_per_step);
	double expected_final_WE = initial_WE - melt_energy_mm_per_dt;
	
	MeltInfo info = firn->melt(m);
	
	water_flux excess = m - info.remaining_energy;

	EXPECT_DOUBLE_EQ(info.melt.value,initial_WE);	

	EXPECT_DOUBLE_EQ(firn->water_equivalent().value,
			0.0);
	
	EXPECT_NEAR(info.remaining_energy.mm_per_dt(p.seconds_per_step),
			m.mm_per_dt(p.seconds_per_step) - initial_WE,MARGIN);
};

TEST_F(LayeredFirnMeltTest,MeltAllLayersExact)
{
	constexpr double initial_WE = std::accumulate(WE.begin(),WE.end(),0.0);
	constexpr auto melt_energy_mm_per_dt = initial_WE;
	water_flux m = water_flux<FluxType::latent>::
		from_mm_per_dt(melt_energy_mm_per_dt,p.seconds_per_step);
	double expected_final_WE = initial_WE - melt_energy_mm_per_dt;
	
	MeltInfo info = firn->melt(m);
	
	water_flux excess = m - info.remaining_energy;

	EXPECT_DOUBLE_EQ(info.melt.value,initial_WE);	

	EXPECT_DOUBLE_EQ(firn->water_equivalent().value,
			0.0);
	
	EXPECT_EQ(info.remaining_energy.mm_per_dt(p.seconds_per_step),
			0.0);

};

/*
 * IceTest: Test fixture for the Ice class, modelling ice as a bulk, single layer.
 */

class IceTest : public ::testing::Test
{
protected:
    Params p = test_default_params();
    std::unique_ptr<Ice> ice;
};

TEST_F(IceTest,ConstructEmpty)
{
    ice = std::make_unique<Ice>(&p);

    EXPECT_EQ(ice->water_equivalent().value,0.0);
};

TEST_F(IceTest,ConstructFromEmptuMilimeterStruct)
{
    Units::Milimeters WE{0.0};
    ice = std::make_unique<Ice>(&p,WE);

    EXPECT_EQ(ice->water_equivalent().value,0.0);
};
TEST_F(IceTest,ConstructsFromMilimetreStruct)
{   
    Units::Milimeters WE{500.34};
    ice = std::make_unique<Ice>(&p,WE);

    EXPECT_EQ(ice->water_equivalent(),WE);
};

TEST_F(IceTest,Accumulate)
{
    Units::Milimeters WE{125.0};
    ice = std::make_unique<Ice>(&p,WE);

    constexpr auto incoming = 340.5;
    WE += Units::Milimeters{incoming};
    ice->accumulate(Units::Milimeters{incoming});

    auto current_WE = ice->water_equivalent();

    EXPECT_EQ(current_WE.value,WE.value);
};

TEST_F(IceTest,MeltSmall)
{
    Units::Milimeters WE{3030.0};
    ice = std::make_unique<Ice>(&p,WE);
    Units::Milimeters to_melt{125.0};

    auto m = water_flux<FluxType::latent>::from_mm_per_dt(to_melt.value,p.seconds_per_step);

    auto melt_info = ice->melt(m);

    EXPECT_EQ(melt_info.melt,to_melt);

    EXPECT_EQ(ice->water_equivalent().value,WE.value - m.mm_per_dt(p.seconds_per_step));

};

TEST_F(IceTest,MeltBig)
{
    Units::Milimeters WE{1245.0};
    ice = std::make_unique<Ice>(&p,WE);
    Units::Milimeters to_melt{WE.value * 1.35};

    auto m = water_flux<FluxType::latent>::from_mm_per_dt(to_melt.value,p.seconds_per_step);

    auto melt_info = ice->melt(m);

    EXPECT_EQ(melt_info.melt,WE);

    EXPECT_DOUBLE_EQ(melt_info.remaining_energy.mm_per_dt(p.seconds_per_step),WE.value*0.35);
};
    
class StateTest : public ::testing::Test
{
protected:
    Params p = test_default_params();

	Ice ice()
	{
		Units::Milimeters ice_WE{100.0};
		return Ice{&p,ice_WE};
	};

	LayeredFirn firn()
	{
		struct Firn_WE
		{
			Layer layer;
			Units::Milimeters mm;
		};

		std::vector<Units::Milimeters> init {{
			{100.0},{200.0},{350.0},{75.0}
		}};

		std::vector<Layer> layers;

		for (auto i : init)
		{
			layers.push_back(i);
		}


		return LayeredFirn{&p,layers};
	}

};

TEST_F(StateTest,ConstructEmpty)
{
    State state(&p);

    auto WE = state.total_water_equiv();
	auto depth = state.total_depth();

    Units::Milimeters expected{0.0};

    EXPECT_DOUBLE_EQ(WE.value,expected.value);
	constexpr auto MM_PER_M = 1000.0;
	EXPECT_DOUBLE_EQ(depth.value * MM_PER_M,expected.value);
};

TEST_F(StateTest,ConstructWithIceOnly)
{
	Ice i = ice();
    State state(&p,i);

	auto IWE = i.water_equivalent();
	auto state_WE = state.total_water_equiv();
	
	constexpr auto M_PER_MM = 1e-3;
	const auto water_density = PhysConst::water_reference_density();
	auto expected_depth = IWE.value * water_density
		/ PhysConst::rho_ice() * M_PER_MM;
	auto depth = state.total_depth();

	EXPECT_DOUBLE_EQ(i.water_equivalent().value,
			state.total_water_equiv().value);
	EXPECT_DOUBLE_EQ(depth.value,expected_depth);
	
};

TEST_F(StateTest,ConstructWithIceAndFirn)
{
	auto i = ice();
	auto f = firn();

	State s(f,i); // No p, not needed by s
				  
	Units::Milimeters init = i.water_equivalent();
	init += f.water_equivalent();

	constexpr auto M_PER_MM = 1e-3;
	const auto water_density = PhysConst::water_reference_density();
	auto expected_depth = i.water_equivalent().value * water_density
		/ PhysConst::rho_ice() * M_PER_MM;
	auto& layers = f.get_layers();
	expected_depth += std::accumulate(layers.begin(),layers.end(),
			0.0,
			[](double acc,Layer layer) {
			acc += layer.height().value;
			return acc;
			});
	auto depth = s.total_depth();


	EXPECT_DOUBLE_EQ(init.value,
			s.total_water_equiv().value);
	EXPECT_DOUBLE_EQ(expected_depth,depth.value);
};

/*
 * DailyMeltTest: Testing of the functions in the Glacier::Melt namespace which contain the 
 * melt logic. Used by the Model class directly.
 */

class DailyMeltTest : public ::testing::Test
{
	std::array<Units::Milimeters,4> WE {{
			{200.0},
			{500.0},
			{250.0},
			{325.0}
		}};
protected:
    bool system_set = false; 
    Params p = test_default_params();
    LayeredFirn firn{&p};
    void set_firn()
    {
		for (auto we : WE)
		{
			firn.accumulate(we);
		}
		
		double sum = std::accumulate(WE.begin(),WE.end(),0.0,
				[](double acc, const Units::Milimeters mm)
				{
					return acc + mm.value;
				});
		EXPECT_NEAR(sum,firn.water_equivalent().value,MARGIN);

		auto& layers = firn.get_layers();

		for (auto layer : layers)
		{
			EXPECT_LT(layer.density().value,p.firn_to_ice_density);
		};

    };
    
    Ice ice{&p};
    void set_ice()
    {
        Units::Milimeters temp{400.0};
        ice.accumulate(temp);
    };

    Units::Milimeters swe{0.0};
    void set_swe()
    {
        swe.value = 450.0;
    }

	water_flux<FluxType::latent> melt_energy = water_flux<FluxType::latent>::from_W_per_m_squared(0.0);
	void SetEnergySmall()
	{	
		double firn_WE = std::accumulate(WE.begin(),WE.end(),0.0,
				[](double acc,Units::Milimeters mm) {
				acc += mm.value; 
				return acc;});
		auto value = water_flux<FluxType::latent>::from_mm_per_dt(firn_WE/2.0,p.seconds_per_step);
		melt_energy += value;
	};

	void SetEnergyBig()
	{
		double firn_WE = std::accumulate(WE.begin(),WE.end(),0.0,
				[](double acc,Units::Milimeters mm) {
				acc += mm.value; 
				return acc;});
		melt_energy = water_flux<FluxType::latent>::from_mm_per_dt(firn_WE*2.0,p.seconds_per_step);
		//melt_energy += value;
	};

    void NoGlacier() {
    };

    void FullGlacier() {
        set_swe();
        set_firn();
        set_ice();
    };

    void no_firn()
    {
        set_swe();
        set_ice();
    };

    void only_snow()
    {
        set_swe();
    };

    void no_snow()
    {
        set_firn();
        set_ice();
    };

    void only_firn()
    {
        set_firn();
    };

    void no_ice()
    {
        set_swe();
        set_firn();
    };
    
    void only_ice()
    {
        set_ice();
    };
	
	template<typename GlacierState>
	void test_scenario(GlacierState f,
			Melt::MeltScenario withoutEnergy,
			Melt::MeltScenario withEnergy)
	{
		f();

		State s(firn,ice);

		Units::Milimeters melt_energy_mm_per_dt{melt_energy.mm_per_dt(p.seconds_per_step)};

		EXPECT_EQ(Melt::get_scenario(s,melt_energy_mm_per_dt,swe), 
				withoutEnergy);

		SetEnergySmall();


		melt_energy_mm_per_dt += Units::Milimeters{melt_energy.mm_per_dt(p.seconds_per_step)};

		EXPECT_EQ(Melt::get_scenario(s,melt_energy_mm_per_dt,swe),
				withEnergy);
	};
	
	template<typename GlacierState>
	void test_scenario(GlacierState f,
			Melt::MeltScenario withoutEnergy)
	{
		test_scenario(f,withoutEnergy,withoutEnergy);
	};
};

TEST_F(DailyMeltTest,FullGlacierIsNoMelt)
{
	test_scenario([this]() { FullGlacier(); },
			Melt::MeltScenario::NoMelt);
};

TEST_F(DailyMeltTest,SWEandFirnIsNoMelt)
{
	test_scenario([this]() { no_ice(); },
			Melt::MeltScenario::NoMelt);
};

TEST_F(DailyMeltTest,OnlySnowIsNoMelt)
{
	test_scenario([this]() { only_snow(); },
			Melt::MeltScenario::NoMelt);
};

TEST_F(DailyMeltTest,NoGlacierIsNoMelt)
{
	test_scenario([this]() { NoGlacier(); },
			Melt::MeltScenario::NoMelt);
};

TEST_F(DailyMeltTest,NoSnowGlacier)
{
	test_scenario([this]() { no_snow(); },
			Melt::MeltScenario::NoMelt,
			Melt::MeltScenario::FirnMelt);

	SetEnergyBig();

	State s(firn,ice);

	Units::Milimeters melt_energy_mm_per_dt{melt_energy.mm_per_dt(p.seconds_per_step)};


	EXPECT_EQ(Melt::get_scenario(s,melt_energy_mm_per_dt,swe), 
			Melt::MeltScenario::FirnAndIceMelt);
};

TEST_F(DailyMeltTest,OnlyIce)
{
	test_scenario([this]() { only_ice(); },
			Melt::MeltScenario::NoMelt,
			Melt::MeltScenario::IceMelt);
};

TEST_F(DailyMeltTest,SnowAndIceNoMelt)
{
	test_scenario([this]() { no_firn(); },
			Melt::MeltScenario::NoMelt);
};

TEST_F(DailyMeltTest,OnlyFirn)
{
	test_scenario([this]() { only_firn(); },
			Melt::MeltScenario::NoMelt,
			Melt::MeltScenario::FirnMelt);
};

TEST_F(DailyMeltTest,NoMelt)
{
	auto scenario = Melt::MeltScenario::NoMelt;

	FullGlacier();

	State s(firn,ice);

	auto melt_energy = water_flux<FluxType::latent>::from_W_per_m_squared(1000.0);
	auto result = Melt::compute_melt(scenario,s,melt_energy);

	EXPECT_FALSE(result);
};

TEST_F(DailyMeltTest,FirnMelt)
{
	auto scenario = Melt::MeltScenario::FirnMelt;

	no_snow();

	State s(firn,ice);
	auto init_FWE = s.firn.water_equivalent();	
	water_flux<FluxType::latent> melt_energy = water_flux<FluxType::latent>::from_mm_per_dt(firn.water_equivalent().value / 2.0,p.seconds_per_step);
	auto result = Melt::compute_melt(scenario,s,melt_energy);
	auto final_FWE = s.firn.water_equivalent();

	EXPECT_DOUBLE_EQ(result->melt.value,melt_energy.mm_per_dt(p.seconds_per_step));
	EXPECT_DOUBLE_EQ(final_FWE.value,init_FWE.value - melt_energy.mm_per_dt(p.seconds_per_step)); 
};

TEST_F(DailyMeltTest,FirnAndIceMelt)
{
	auto scenario = Melt::MeltScenario::FirnAndIceMelt;

	no_snow();

	State s(firn,ice);
	auto init_IWE = s.ice.water_equivalent();
	constexpr auto ice_melt_fraction = 0.5;
	
	water_flux<FluxType::latent> melt_energy = water_flux<FluxType::latent>::from_mm_per_dt(firn.water_equivalent().value +init_IWE.value * ice_melt_fraction,p.seconds_per_step);
	auto result = Melt::compute_melt(scenario,s,melt_energy);
	auto final_FWE = s.firn.water_equivalent();
	auto& layers = s.firn.get_layers();
	auto final_IWE = s.ice.water_equivalent();

	EXPECT_DOUBLE_EQ(result->melt.value,melt_energy.mm_per_dt(p.seconds_per_step));
	EXPECT_DOUBLE_EQ(final_FWE.value,0.0);
	EXPECT_EQ(layers.size(),0);
	EXPECT_DOUBLE_EQ(final_IWE.value,init_IWE.value * ice_melt_fraction);

};

TEST_F(DailyMeltTest,FirnAndIceMeltAll)
{
	auto scenario = Melt::MeltScenario::FirnAndIceMelt;

	no_snow();

	State s(firn,ice);
	auto init_IWE = s.ice.water_equivalent();
	constexpr auto ice_melt_fraction = 1.5;
	
	water_flux<FluxType::latent> melt_energy = water_flux<FluxType::latent>::from_mm_per_dt(firn.water_equivalent().value +init_IWE.value * ice_melt_fraction,p.seconds_per_step);
	auto result = Melt::compute_melt(scenario,s,melt_energy);
	auto final_FWE = s.firn.water_equivalent();
	auto& layers = s.firn.get_layers();
	auto final_IWE = s.ice.water_equivalent();

	EXPECT_DOUBLE_EQ(result->melt.value,melt_energy.mm_per_dt(p.seconds_per_step) - init_IWE.value * 0.5);
	EXPECT_DOUBLE_EQ(final_FWE.value,0.0);
	EXPECT_EQ(layers.size(),0);
	EXPECT_DOUBLE_EQ(final_IWE.value,0.0);
};

TEST_F(DailyMeltTest,IceMelt)
{
	auto scenario = Melt::MeltScenario::IceMelt;

	only_ice();

	State s(firn,ice);

	auto init_IWE = s.ice.water_equivalent();

	water_flux<FluxType::latent> melt_energy = 
		water_flux<FluxType::latent>::from_mm_per_dt(s.ice.water_equivalent().value * 0.75,p.seconds_per_step);
	auto result = Melt::compute_melt(scenario,s,melt_energy);

	auto final_IWE = s.ice.water_equivalent();

	EXPECT_EQ(result->melt.value,melt_energy.mm_per_dt(p.seconds_per_step));
	EXPECT_EQ(final_IWE.value,init_IWE.value - melt_energy.mm_per_dt(p.seconds_per_step));

};
/*
 * The Updater is responsible for the infrequent updating part 
 * of the glacier season. The following processes occur:
 *
 * 1. Converts SWE to firn.
 * 2. Converts firn to ice.
 * 3. Increases the density of firn.
 *
 * The following is a test for the Updater class, constructed only
 * as needed on the stack, obtaining references to firn and ice 
 * objects, as well as the swe as Units::Milimeter instance.
 *
 * The best approach is to have a member function for each:
 *
 * 1. SWE_to_firn()
 * 2. firn_to_ice()
 * 3. densify_firn()
 *
 */
class YearlyUpdaterTest : public ::testing::Test
{
protected:
	Params p = test_default_params();
	Ice ice;
	
	YearlyUpdaterTest() : ice(&p) {};
	
	// Construct fresh firn, no layers high density
	LayeredFirn fresh_firn() {
		std::vector<Units::Milimeters> WE {{
			{200.0},
			{500.0},
			{250.0},
			{325.0}
		}};

		LayeredFirn firn(&p);

		for (auto we : WE)
		{
			firn.accumulate(we);
		}
		
		double sum = std::accumulate(WE.begin(),WE.end(),0.0,
				[](double acc, const Units::Milimeters mm)
				{
					return acc + mm.value;
				});
		EXPECT_NEAR(sum,firn.water_equivalent().value,MARGIN);

		auto& layers = firn.get_layers();

		for (auto layer : layers)
		{
			EXPECT_LT(layer.density().value,p.firn_to_ice_density);
		};

		return firn;
	};

	//construct old firn, several (but not all) layers with 
	//high density
	LayeredFirn old_firn()
	{
		struct Local
		{
			Layer::Height height;
			Layer::Density density;
		};
		std::array<Local,6> local{{
			{{4.0},{250.0}},{{3.2},{420.0}},{{2.5},{750.0}},
			{{1.2},{p.firn_to_ice_density + 20.0}},
			{{0.6},{p.firn_to_ice_density + 45.0}},
			{{0.3},{p.firn_to_ice_density + 100.0}}
		}};
			
		std::deque<Layer> layers;

		for (auto l : local)
		{
			Layer layer(l.height,l.density);
			layers.push_front(layer);
		}

		LayeredFirn firn_local(&p,layers);

		auto& layers_ref = firn_local.get_layers();
		EXPECT_EQ(layers_ref.size(),6);
		
		size_t count_over_critical = 0;
		for (auto layer : layers_ref)
		{
			if (layer.density().value > p.firn_to_ice_density)
				++count_over_critical;
		}
		EXPECT_EQ(count_over_critical,3);
		
		return firn_local;
	};
		
	void some_ice()
	{
		Units::Milimeters ice_init_WE{6e3};
		ice.accumulate(ice_init_WE);
	};

};

TEST_F(YearlyUpdaterTest,NoSWEFreshFirnUnchanged)
{
	LayeredFirn firn = fresh_firn();
	Units::Milimeters swe{0.0};
    auto swe_copy = swe;
	auto layers_init = firn.get_layers();
	
	Updates::swe_to_firn(p,firn,swe);

	auto layers_after = firn.get_layers();

	EXPECT_EQ(layers_init.size(),layers_after.size());
	EXPECT_EQ(layers_init.back().water_equivalent()
            ,layers_after.back().water_equivalent());
    EXPECT_EQ(swe_copy.value,swe.value);
};

TEST_F(YearlyUpdaterTest,SWEupdatesFirn)
{
	LayeredFirn firn = fresh_firn();
	Units::Milimeters swe{150.0};
    auto swe_copy = swe;
	const size_t init_firn_layers =
		firn.get_layers().size();

	Updates::swe_to_firn(p,firn,swe);

	auto& layers = firn.get_layers();

	EXPECT_EQ(layers.size(),init_firn_layers+1);

	EXPECT_NEAR(layers.back().water_equivalent().value,
			swe_copy.value,MARGIN);

    EXPECT_EQ(swe.value,0.0);
};

TEST_F(YearlyUpdaterTest,FreshFirnNotConvertedToIce)
{
	LayeredFirn firn = fresh_firn();
	Units::Milimeters swe{0.0};

    //LayeredFirn firn2{Params()};
	const size_t init_firn_layers =
		firn.get_layers().size();
	const Units::Milimeters init_ice = ice.water_equivalent();

	Updates::firn_to_ice(p,firn,ice);

	EXPECT_EQ(firn.get_layers().size(),init_firn_layers);

	EXPECT_EQ(Units::Milimeters{0.0},ice.water_equivalent());

};

TEST_F(YearlyUpdaterTest,OldFirnConvertsFrontToIce)
{
	LayeredFirn firn = old_firn();
	some_ice();
	Units::Milimeters swe{0.0};

	auto init_layers = firn.get_layers();
	const size_t init_firn_layers =
		init_layers.size();
	const Units::Milimeters init_ice = ice.water_equivalent();

	Updates::firn_to_ice(p,firn,ice);

	auto& new_layers = firn.get_layers();	

	auto it_new = new_layers.rbegin();
	auto it_init = init_layers.rbegin();

	while (it_new != new_layers.rend())
	{
		EXPECT_EQ(it_new->water_equivalent().value,
		it_init->water_equivalent().value);
		++it_new;
		++it_init;
	}
	EXPECT_TRUE(it_init == init_layers.rend()-1);

	EXPECT_EQ(init_layers.back().water_equivalent(),
            new_layers.back().water_equivalent());

	EXPECT_EQ(init_firn_layers-1,new_layers.size());

	EXPECT_NE(init_layers.front().water_equivalent(),
            new_layers.front().water_equivalent());

	EXPECT_EQ(new_layers.front().water_equivalent(),
            init_layers[1].water_equivalent());

	EXPECT_EQ(ice.water_equivalent().value,
			init_ice.value + 
			init_layers.front().water_equivalent().value);
};

TEST_F(YearlyUpdaterTest,FreshFirnUnchangedIceAndFirn)
{
	LayeredFirn firn = fresh_firn();
	some_ice();

	LayeredFirn init_firn = firn;
	const Units::Milimeters init_ice = ice.water_equivalent();

	Updates::firn_to_ice(p,firn,ice);

	EXPECT_EQ(init_firn.water_equivalent().value,
			firn.water_equivalent().value);

	EXPECT_EQ(init_ice.value,ice.water_equivalent().value);
};

TEST_F(YearlyUpdaterTest,ManyIceUpdates)
{
	LayeredFirn firn = old_firn();
	const Units::Milimeters init_ice = ice.water_equivalent();

	auto& layers = firn.get_layers();
    size_t init_layer_num = layers.size();
	auto init_layers = layers;
	Units::Milimeters sum{init_ice};
    size_t iter_tracker = 0;
	while (iter_tracker < init_layer_num)
	{	
		if (layers.front().water_equivalent().value < p.firn_to_ice_density)
			break;
        ++iter_tracker;
		sum += layers.front().water_equivalent();
		Updates::firn_to_ice(p,firn,ice);
		EXPECT_EQ(sum.value,ice.water_equivalent().value);

	}

};

class TestData
{
    State _state;

    // Inputs (set before calling execute_impl)
    Units::Milimeters _swe{0.0};
    Units::Watts_per_m2 _melt_energy{0.0};
    Units::Celsius    _glacier_temperature{273.15};
    bool              _update_now{false};

    // Outputs (recorded by execute_impl)
    double _glacier_water_equivalent{0.0};
    double _total_depth{0.0};
    double _firn_melt{0.0};
    double _ice_melt{0.0};

public:
    explicit TestData(const Params* p) : _state(p) { check_params(p);}
    TestData(const Params* p, const LayeredFirn& f) 
        : _state(p, f) {check_params(p);}  // State(Params*, LayeredFirn&)
    TestData(const Params* p, const Ice& i) 
        : _state(p, i) {check_params(p);}
    TestData(const LayeredFirn& f, const Ice& i) 
        : _state(f, i) {}

	void check_params(const Params* p)
	{
		EXPECT_EQ(p->seconds_per_step,3600) << "See hard-coded value in set_melt_energy";
	};

    // ── Setters for test configuration ──
    void set_swe(double v) { _swe.value = v; }
    void set_melt_energy(double v) { 
		auto wf = water_flux<FluxType::latent>::from_mm_per_dt(v,3600);
		_melt_energy.value = wf.W_per_m_squared(); }
    void set_glacier_temperature(double v) { _glacier_temperature.value = v; }
    void set_update_now(bool v) { _update_now = v; }

    // ── GlacierData concept interface ──
    Units::Milimeters swe() { return _swe; }
    Units::Watts_per_m2 melt_energy() { return _melt_energy; }
    Units::Celsius glacier_temperature() { return _glacier_temperature; }
    bool update_now() { return _update_now; }

    void glacier_water_equivalent(double v) { _glacier_water_equivalent = v; }
    void total_depth(double v) { _total_depth = v; }
    void firn_melt(double v) { _firn_melt = v; }
    void ice_melt(double v) { _ice_melt = v; }

    // ── Extra methods used by execute_impl ──
    State& get_state() { return _state; }
    const State& get_state() const { return _state; }
    Units::Milimeters snowmelt() { return Units::Milimeters{0.0}; } // unused by logic
    Units::Milimeters rainfall() { return Units::Milimeters{0.0}; } // unused by logic

    // ── Getters for test assertions ──
    double recorded_glacier_water_equivalent() const { return _glacier_water_equivalent; }
    double recorded_total_depth() const { return _total_depth; }
    double recorded_firn_melt() const { return _firn_melt; }
    double recorded_ice_melt() const { return _ice_melt; }
};

static_assert(GlacierData<TestData>);

class ModelTest : public ::testing::Test
{
protected:
    Params p;
    Model<TestData> model;

    // Helper: run model on data
    void run(TestData& d) { model.execute_impl(d); }

    // Helper: create firn with known WE
    LayeredFirn make_firn(double we_mm)
    {
        LayeredFirn f(&p);
        if (we_mm > 0.0)
            f.accumulate(Units::Milimeters{we_mm});
        return f;
    }

    // Helper: create ice with known WE
    Ice make_ice(double we_mm)
    {
        return Ice(&p, Units::Milimeters{we_mm});
    }
};

TEST_F(ModelTest, SwePositive_EmptyGlacier_NoMelt)
{
    TestData d(&p);
    d.set_swe(100.0);
    d.set_melt_energy(500.0);  // would melt if SWE were 0

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, SwePositive_WithFirn_NoMelt)
{
    auto firn = make_firn(200.0);
    TestData d(firn, make_ice(0.0));
    d.set_swe(50.0);
    d.set_melt_energy(1000.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, SwePositive_WithIce_NoMelt)
{
    TestData d(&p, make_ice(500.0));
    d.set_swe(10.0);
    d.set_melt_energy(200.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, SwePositive_WithFirnAndIce_NoMelt)
{
    auto firn = make_firn(300.0);
    TestData d(firn, make_ice(400.0));
    d.set_swe(1.0);
    d.set_melt_energy(9999.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, ZeroMeltEnergy_EmptyGlacier_NoMelt)
{
    TestData d(&p);
    d.set_swe(0.0);
    d.set_melt_energy(0.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, ZeroMeltEnergy_WithFirn_NoMelt)
{
    auto firn = make_firn(200.0);
    TestData d(firn, make_ice(0.0));
    d.set_swe(0.0);
    d.set_melt_energy(0.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, ZeroMeltEnergy_WithIce_NoMelt)
{
    TestData d(&p, make_ice(500.0));
    d.set_swe(0.0);
    d.set_melt_energy(0.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, FirnOnlyMelt_PartialMelt)
{
    // Firn only, no ice → FirnMelt scenario
    auto firn = make_firn(500.0);
    TestData d(firn, make_ice(0.0));
    d.set_swe(0.0);
    d.set_melt_energy(100.0);  // less than firn WE
    d.set_glacier_temperature(0.0);

    run(d);

    EXPECT_GT(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
    // Firn WE should decrease
    EXPECT_LT(d.get_state().firn.water_equivalent().value, 500.0);
}

TEST_F(ModelTest, IceOnlyMelt)
{
    // Ice only, no firn → IceMelt scenario
    TestData d(&p, make_ice(500.0));
    d.set_swe(0.0);
    d.set_melt_energy(100.0);
    d.set_glacier_temperature(0.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    EXPECT_GT(d.recorded_ice_melt(), 0.0);
    EXPECT_LT(d.get_state().ice.water_equivalent().value, 500.0);
}

TEST_F(ModelTest, FirnAndIceMelt_MeltExceedsFirn)
{
    // Both firn and ice present, melt energy exceeds firn WE
    // → FirnAndIceMelt scenario
    auto firn = make_firn(50.0);  // small firn
    auto ice = make_ice(500.0);
    TestData d(firn, ice);
    d.set_swe(0.0);
    d.set_melt_energy(200.0);  // well exceeds firn
    d.set_glacier_temperature(0.0);

    double initial_firn_we = d.get_state().firn.water_equivalent().value;
    double initial_ice_we = d.get_state().ice.water_equivalent().value;

    run(d);

    // In FirnAndIceMelt, firn_melt should equal initial firn WE (all firn melted)
    EXPECT_NEAR(d.recorded_firn_melt(), initial_firn_we, MARGIN);
    // Ice melt should be positive
    EXPECT_GT(d.recorded_ice_melt(), 0.0);
    // Ice melt = initial_ice_we - remaining ice
    double expected_ice_melt = initial_ice_we - d.get_state().ice.water_equivalent().value;
    EXPECT_NEAR(d.recorded_ice_melt(), expected_ice_melt, MARGIN);
}

TEST_F(ModelTest, FirnAndIce_MeltDoesNotExceedFirn)
{
    // Both present but melt energy < firn WE → FirnMelt only
    auto firn = make_firn(500.0);
    auto ice = make_ice(500.0);
    TestData d(firn, ice);
    d.set_swe(0.0);
    d.set_melt_energy(50.0);  // less than firn
    d.set_glacier_temperature(0.0);

    run(d);

    EXPECT_GT(d.recorded_firn_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_ice_melt(), 0.0);
}

TEST_F(ModelTest, UpdateNowTrue_SweConvertedToFirn)
{
    TestData d(&p);
    d.set_swe(300.0);
    d.set_melt_energy(0.0);
    d.set_update_now(true);

    EXPECT_EQ(d.get_state().firn.get_layers().size(), 0u);

    run(d);

    // SWE should have been accumulated as a new firn layer
    EXPECT_GT(d.get_state().firn.get_layers().size(), 0u);
    EXPECT_GT(d.get_state().firn.water_equivalent().value, 0.0);
}

TEST_F(ModelTest, UpdateNowFalse_SweNotConvertedToFirn)
{
    TestData d(&p);
    d.set_swe(300.0);
    d.set_melt_energy(0.0);
    d.set_update_now(false);

    run(d);

    // No firn conversion should happen
    EXPECT_EQ(d.get_state().firn.get_layers().size(), 0u);
}

TEST_F(ModelTest, UpdateNowTrue_DenseFirnConvertsToIce)
{
    // Manually create a firn layer with density above firn_to_ice threshold
    Layer::Height h{1.0};
    Layer::Density rho{p.firn_to_ice_density + 10.0};  // above threshold
    Layer dense_layer(h, rho);

    LayeredFirn firn(&p, std::vector<Layer>{dense_layer});
    TestData d(firn, Ice(&p));
    d.set_swe(0.0);
    d.set_melt_energy(0.0);
    d.set_update_now(true);

    double initial_firn_we = d.get_state().firn.water_equivalent().value;
    EXPECT_GT(initial_firn_we, 0.0);
    EXPECT_DOUBLE_EQ(d.get_state().ice.water_equivalent().value, 0.0);

    run(d);

    // Dense layer should have been converted to ice
    EXPECT_GT(d.get_state().ice.water_equivalent().value, 0.0);
}

TEST_F(ModelTest, UpdateNowFalse_DenseFirnDoesNotConvertToIce)
{
    Layer::Height h{1.0};
    Layer::Density rho{p.firn_to_ice_density + 10.0};
    Layer dense_layer(h, rho);

    LayeredFirn firn(&p, std::vector<Layer>{dense_layer});
    TestData d(firn, Ice(&p));
    d.set_swe(0.0);
    d.set_melt_energy(0.0);
    d.set_update_now(false);

    run(d);

    // No conversion
    EXPECT_DOUBLE_EQ(d.get_state().ice.water_equivalent().value, 0.0);
    EXPECT_GT(d.get_state().firn.water_equivalent().value, 0.0);
}

TEST_F(ModelTest, OutputsMatchStateAfterExecution_EmptyGlacier)
{
    TestData d(&p);
    d.set_swe(0.0);
    d.set_melt_energy(0.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_glacier_water_equivalent(), 
                     d.get_state().total_water_equiv().value);
    EXPECT_DOUBLE_EQ(d.recorded_total_depth(),
                     d.get_state().total_depth().value);
}

TEST_F(ModelTest, OutputsMatchStateAfterExecution_WithMelt)
{
    auto firn = make_firn(300.0);
    TestData d(firn, make_ice(200.0));
    d.set_swe(0.0);
    d.set_melt_energy(50.0);
    d.set_glacier_temperature(0.0);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_glacier_water_equivalent(),
                     d.get_state().total_water_equiv().value);
    EXPECT_DOUBLE_EQ(d.recorded_total_depth(),
                     d.get_state().total_depth().value);
}

TEST_F(ModelTest, OutputsMatchStateAfterExecution_WithUpdate)
{
    TestData d(&p);
    d.set_swe(500.0);
    d.set_melt_energy(0.0);
    d.set_update_now(true);

    run(d);

    EXPECT_DOUBLE_EQ(d.recorded_glacier_water_equivalent(),
                     d.get_state().total_water_equiv().value);
    EXPECT_DOUBLE_EQ(d.recorded_total_depth(),
                     d.get_state().total_depth().value);
}

TEST_F(ModelTest, DifferentTemperatures_DifferentMeltAmounts)
{
    // Same melt_energy but different temperatures should produce 
    // different water_flux conversions → different melt amounts
    auto firn1 = make_firn(500.0);
    TestData d1(firn1, make_ice(0.0));
    d1.set_swe(0.0);
    d1.set_melt_energy(100.0);
    d1.set_glacier_temperature(-10.0);

    auto firn2 = make_firn(500.0);
    TestData d2(firn2, make_ice(0.0));
    d2.set_swe(0.0);
    d2.set_melt_energy(100.0);
    d2.set_glacier_temperature(-30.0);

    run(d1);
    run(d2);

    // The melt amounts should differ due to temperature-dependent conversion
    // (If water_flux::from_W_per_m_squared uses temperature, they'll differ.
    //  If not, this test documents that temperature has no effect on melt amount.)
    // Either way, both should report zero ice melt
    EXPECT_DOUBLE_EQ(d1.recorded_ice_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d2.recorded_ice_melt(), 0.0);

	EXPECT_NE(d1.recorded_firn_melt(),d2.recorded_firn_melt());
}

TEST_F(ModelTest, MeltExceedsTotalGlacier_FullDepletion)
{
    auto firn = make_firn(50.0);
    auto ice = make_ice(50.0);
    TestData d(firn, ice);
    d.set_swe(0.0);
    d.set_melt_energy(99999.0);  // massive melt
    d.set_glacier_temperature(0.0);

    run(d);

    // Everything should be melted
    EXPECT_NEAR(d.get_state().total_water_equiv().value, 0.0, 1e-6);
    EXPECT_DOUBLE_EQ(d.recorded_glacier_water_equivalent(), 
                     d.get_state().total_water_equiv().value);
}

TEST_F(ModelTest, IceMeltThenUpdate_NoFirnConversionOnEmptyFirn)
{
    // Ice-only glacier, melt some ice, then update
    TestData d(&p, make_ice(500.0));
    d.set_swe(0.0);
    d.set_melt_energy(50.0);
    d.set_glacier_temperature(0.0);
    d.set_update_now(true);

    run(d);

    EXPECT_GT(d.recorded_ice_melt(), 0.0);
    EXPECT_DOUBLE_EQ(d.recorded_firn_melt(), 0.0);
    // update_now with no SWE and no dense firn → no state change from updates
    EXPECT_EQ(d.get_state().firn.get_layers().size(), 0u);
}
