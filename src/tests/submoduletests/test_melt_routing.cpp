#include "gtest/gtest.h"
#include <cstdlib>
#include <stdexcept>
#include "melt_routing_glacier.hpp"
#include "PhysConst.h"

namespace Units = PhysConst::units;
using namespace GlacierRouting;

inline Params test_params()
{
    Params p;

    p.chain_length = 10u;
    p.k = 4e3;
    p.seconds_per_step = 3600u;

	return p;
};

class LinearReservoirTest : public ::testing::Test
{
protected:
	void run(const double k, const size_t dt, std::vector<double> inputs) {
		LinearReservoir res(k,dt);
		auto c0 = dt / (2*k + dt);
		auto c1 = (2*k - dt) / (2*k + dt);

		auto last_input = 0.0;
		auto last_output = 0.0;

		for (auto input : inputs)
		{
			auto result = res.step(input);
			auto expected_result = c0 * (input + last_input) 
				+ c1 * last_output;
			EXPECT_DOUBLE_EQ(result,expected_result);

			last_input = input;
			last_output = expected_result;
		}
	}
};

TEST_F(LinearReservoirTest,ZeroKInIsOut)
{
	LinearReservoir res(0.0,1.0);

	std::vector<double> Vec{105.0,250.0,1234.111,24591.0};
	for (auto v : Vec)
	{
		auto output = res.step(v);
		EXPECT_EQ(output,v);
	}
};

TEST_F(LinearReservoirTest,NonZeroKLessThanCriticalThrows)
{
	const auto dt = 1000.0;
	const auto k = dt / 2 - 0.01;
	auto func = [=]() {
		LinearReservoir res(k,dt);
	};

	EXPECT_THROW(func(),std::runtime_error);
};

TEST_F(LinearReservoirTest,StepTest)
{
	const auto dt = 1500.0;
	const auto k = 3050.0;

	std::vector<double> inputs{1.0,2.5,1.3,8.2};

	run(k,dt,inputs);
};

TEST_F(LinearReservoirTest,StepTest2)
{
	size_t N = 1000;
	std::vector<double> inputs;
	auto dt = 3600u;
	auto k = 5000.0;
	auto amplitude = 250.5;
	for (size_t i = 0; i < N; ++i)
	{
		auto val = amplitude * 
			( std::sin( amplitude* (i * dt)) + 1 );

		inputs.emplace_back(val);
	}

	run(k,dt,inputs);

}

class GlacierReservoirTest : public ::testing::Test
{
protected:
	Params p = test_params();

};

TEST_F(GlacierReservoirTest,ZeroUntilChainLengthSteps)
{
	GlacierReservoir res(&p);
	
	double initial_input = 150.0;

	res.step(initial_input);

	for (size_t i = 1; i < p.chain_length -1; ++i)
	{
		auto result = res.step(0.0);	
		EXPECT_EQ(result,0.0);
	}

	auto result = res.step(0.0);
	EXPECT_GT(result,0.0);
};

TEST_F(GlacierReservoirTest,ChainLengthOneSameAsLinearReservoir)
{
	p.chain_length = 1;
	LinearReservoir lin_res(p.k,p.seconds_per_step);
	GlacierReservoir gla_res(&p);

	double input = 100.0;

	auto out_lin = lin_res.step(input);
	auto out_gla = gla_res.step(input);

	EXPECT_EQ(out_lin,out_gla);

	// Since they are equal, this just makes sure they aren't both
	// zero
	EXPECT_GT(out_lin + out_gla,0.0);
};

TEST_F(GlacierReservoirTest,LongChainSameAsLinearReservoir)
{
	p.chain_length = 75u;
	LinearReservoir lin_res(p.k,p.seconds_per_step);
	GlacierReservoir gla_res(&p);

	size_t N = 500;
	std::vector<double> inputs;
	auto dt = p.seconds_per_step;
	auto amplitude = 250.5;
	for (size_t i = 0; i < N; ++i)
	{
		auto val = amplitude * 
			( std::sin( amplitude* (i * dt)) + 1 );

		inputs.emplace_back(val);
	}

	for (size_t i = 0; i < N; ++i)
	{
		auto out_lin = 0.0;

		if ( i + 1 >= p.chain_length)	
			out_lin = lin_res.step(inputs[i - p.chain_length + 1]);

		auto out_gla = gla_res.step(inputs[i]);
		
		if (i < p.chain_length - 1)
		{
			EXPECT_DOUBLE_EQ(out_gla,0.0);
		}
		else
			ASSERT_DOUBLE_EQ(out_gla,out_lin);
	};
};

TEST(GlacierRoutingStateTest,Construct)
{
	Params p = test_params();

	State s(&p);
};

class Data
{
	State s;
	struct Inout
	{
		double melt = 0.0;
		double delayed = 0.0;
	};
	Inout snow;
	Inout firn;
	Inout ice;
	double _total_delayed;
public:
    State& get_state() { return s; }
	Data(const Params* _p) : s(_p) {};
    const Units::Milimeters snowmelt() 
	{ return Units::Milimeters{snow.melt}; }
    const Units::Milimeters firnmelt()
	{ return Units::Milimeters{firn.melt}; }
    const Units::Milimeters icemelt()
	{ return Units::Milimeters{ice.melt}; }

	void snowmelt(const double T) { snow.melt = T; };
	void firnmelt(const double T) { firn.melt = T; };
	void icemelt(const double T) { ice.melt = T; };


    void snowmelt_delayed(const double T) 
	{ snow.delayed = T; };
    void firnmelt_delayed(const double T) 
	{ firn.delayed = T; };
    void icemelt_delayed(const double T) 
	{ ice.delayed = T; };
    void total_delayed(const double T) 
	{ _total_delayed = T; };

	const double snowmelt_delayed()
	{ return snow.delayed; }
	const double firnmelt_delayed()
	{ return firn.delayed; }
	const double icemelt_delayed()
	{ return ice.delayed; }
};

class GlacierRoutingModelTest : public ::testing::Test
{
protected:
	Model<Data> model;
	Params& p = model.get_params();
	std::vector<double> snowmelts;
	std::vector<double> firnmelts;
	std::vector<double> icemelts;
	const size_t N = 500u;
	void set_up_inputs(const double amplitude,std::vector<double>& melts)
	{
		auto dt = p.seconds_per_step;
		for (size_t i = 0; i < N; ++i)
		{
			auto val = amplitude * 
				( std::sin( amplitude* (i * dt)) + 1 );

			melts.emplace_back(val);
		}
	};

	void SetUp() override {
		p = test_params();
		p.chain_length = 10u;
	};
};

TEST_F(GlacierRoutingModelTest,ExecuteImpl)
{
	set_up_inputs(150.0,snowmelts);
	set_up_inputs(45.3,firnmelts);
	set_up_inputs(12.3,icemelts);
	Data d(&p);

	struct Res
	{
		LinearReservoir snow;
		LinearReservoir firn;
		LinearReservoir ice;
		Res(const double k, const size_t dt)
			: snow(k,dt), firn(k,dt), ice(k,dt) {};
	} res(p.k,p.seconds_per_step);

	for (size_t i = 0; i < N; ++i)
	{
		d.snowmelt(snowmelts[i]);
		d.firnmelt(firnmelts[i]);
		d.icemelt(icemelts[i]);

		model.execute(d);

		auto out_snow = 0.0;
		auto out_firn = 0.0;
		auto out_ice = 0.0;

		if ( i + 1 >= p.chain_length)	
		{
			out_snow = res.snow.step(snowmelts[i - p.chain_length + 1]);
			out_firn = res.firn.step(firnmelts[i - p.chain_length + 1]);
			out_ice = res.ice.step(icemelts[i - p.chain_length + 1]);
		}

		
		if (i < p.chain_length - 1)
		{
			EXPECT_DOUBLE_EQ(d.snowmelt_delayed(),0.0);
			EXPECT_DOUBLE_EQ(d.firnmelt_delayed(),0.0);
			EXPECT_DOUBLE_EQ(d.icemelt_delayed(),0.0);
		}
		else
		{
			EXPECT_DOUBLE_EQ(d.snowmelt_delayed(),out_snow);
			EXPECT_DOUBLE_EQ(d.firnmelt_delayed(),out_firn);
			EXPECT_NE(d.snowmelt_delayed(),d.firnmelt_delayed());
			ASSERT_DOUBLE_EQ(d.icemelt_delayed(),out_ice);
		}

	}
};




