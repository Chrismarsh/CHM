#include "Ayers.hpp"
#include "gtest/gtest.h"
#include "Soil.h"
/*
 * CrackTest: Wrapper class for tests
 * CrackTest is effectively a mock of Infil_All module but done indirectly. Due to the complexity of the module classes, it was easier to write this.  
 * The member variables with the _ prefix are inputs that are supplied to the constructor of Crack.
 * Default values are given and used for most tests.
 * Other member variables are parameters that are also supplied to Crack unless it has the const specifier, then it is just useful for these tests.
 * Member functions are just tools to enable the tests.
 * Initialization of CrackTest assumes that the frozen period has just begun. 
 * 
 */
class AyersTest : public testing::Test
{
protected:

    AyersTest()
    {
    };

	typedef Soil::soils_na S;

	typedef Ayers<S,&S::ayers_texture> MyAyers;

    double _snowmelt = 1.0;
    double _rainfall = 0.0;
	std::string _ground_cover = "bare_soil";
	std::string _texture = "coarse_over_coarse";		

	S& soils = Soil::get_soil_obj<S>();

	MyAyers DoAyers(double& snowmelt, double& rainfall)
	{
		MyAyers ayers(rainfall,snowmelt,_texture,_ground_cover, soils);

        ayers.run();
		return ayers;
	};

	void DoAssert(MyAyers& ayers,const double inf,const double runoff,const double snowinf)
	{
		ASSERT_EQ(ayers.get_inf(),inf);
		ASSERT_EQ(ayers.get_runoff(),runoff);
		ASSERT_EQ(ayers.get_snow_inf(),snowinf);
	};
};

TEST_F(AyersTest, ZeroInputs)
{
	_rainfall = 0.0;
	_snowmelt = 0.0;

	MyAyers ayers = DoAyers(_snowmelt, _rainfall);

	DoAssert(ayers,0.0,0.0,0.0);
};

TEST_F(AyersTest, NonZeroRainfallSmall)
{
	_rainfall = 1e-6;
	_snowmelt = 0.0;

	MyAyers ayers = DoAyers(_snowmelt,_rainfall);

	DoAssert(ayers,_rainfall,0.0,0.0);

};

TEST_F(AyersTest, NonZeroBigRainfall)
{
	_rainfall = 1e4;
	_snowmelt = 0.0;

	MyAyers ayers = DoAyers(_snowmelt,_rainfall);

	DoAssert(ayers,7.6,_rainfall - 7.6,0.0);

};

TEST_F(AyersTest, NonZeroSnowMelt)
{
	_rainfall = 0.0;
	_snowmelt = 1.3e2;

	MyAyers ayers1 = DoAyers(_snowmelt,_rainfall);

	DoAssert(ayers1,_snowmelt,0.0,_snowmelt);

	_snowmelt = 3.14159;

	MyAyers ayers2 = DoAyers(_snowmelt,_rainfall);

	DoAssert(ayers2,_snowmelt,0.0,_snowmelt);

};

TEST_F(AyersTest, MeltAndRain)
{
	_rainfall = 99;
	_snowmelt = 13;
	
	_ground_cover = "small_grains";
	_texture = "medium_over_medium";

	MyAyers ayers = DoAyers(_snowmelt,_rainfall);
	double maxinfil = 10.2;
	DoAssert(ayers,maxinfil+_snowmelt,_rainfall - maxinfil,_snowmelt);

};


