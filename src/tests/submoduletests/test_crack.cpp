#include "Crack.hpp"
#include "gtest/gtest.h"
#include "CSVreader.hpp"
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
class CrackTest : public testing::Test
{
protected:

    CrackTest()
    {
        status.init();
        status.begin_freeze();
        set_steps_per_day(seconds_per_hour);
    };

    Crack::info status;

    double _snowmelt = 1.0;
    double _rainfall = 0.0;
    double _swe = 100.0;
    double _soil_storage_at_freeze = 50.0;
    double _airtemp = 10.0;
    bool _newday = false;

    double major = 5;
    double min_swe_to_freeze = 25;
    unsigned int infDays = 6;
    bool AllowPriorInf = true;
    double lenstemp = -10.0;
    static constexpr double seconds_per_hour = 3600.0;
    double steps_per_day;

    const double diff = 1e-5;
    
    void set_steps_per_day(double seconds_per_step)
    {
        steps_per_day = 86400 / seconds_per_step;
    };

    void set_newday(bool newday)
    {
        _newday = newday;
    };
    
    // Crack is designed to be declared, initialized and run on every timestep. This function does a single step, while taking a snowmelt as an input.
    Crack run_a_step(double snowmelt)
    {
        Crack model(major,min_swe_to_freeze,infDays,AllowPriorInf,lenstemp,steps_per_day,status);
        
        _swe -= snowmelt;
        model.init_inputs(snowmelt,_rainfall,_swe,_soil_storage_at_freeze,_airtemp,_newday);

        model.run();
        return model;
    };

    // similar to the above but allows the `soil_storage_at_freeze` to be set. 
    Crack restricted_or_unlimited(double melt, double storage)
    {
        const double day_melt = melt; 
        set_newday(true);
        _soil_storage_at_freeze = storage;    
        EXPECT_TRUE(storage == 100.0 || storage == 0.0);
        status.daily_melt_total = day_melt;

        return run_a_step(_snowmelt);
    };

    // Expected infiltration for comparison.
    double expected_inf(double melt)
    {
        return std::min(melt * status.index,status.index / infDays * status.init_SWE);
    };

};


TEST_F(CrackTest, AccumulatesInputsUntilNewDay) {
    
    for (int i = 0; i < steps_per_day; i++)
    {
        Crack model = run_a_step(_snowmelt);

        // inf should be 0.0 because the snow-covered period 
        // just stated and the Crack model applies melt from 
        // day i as infiltration on day i+1
        EXPECT_EQ(model.get_inf(),0.0);
        // checking that snowmelt is accumulating
        EXPECT_EQ(status.daily_melt_total, _snowmelt * (i+1));
    }

    set_newday(true);
    Crack model = run_a_step(_snowmelt);
    // Checks that it restarts properly on a new day
    EXPECT_EQ(status.daily_melt_total, _snowmelt);
    
    
};
TEST_F(CrackTest, ComputeInfOnDailyTotal) {
    set_newday(true);
    double yesterday_melt = _snowmelt * steps_per_day / 2.0;
    status.daily_melt_total = yesterday_melt;
    Crack model = run_a_step(_snowmelt);
    
    // Check that the initial SWE is recorded correctly
    EXPECT_EQ(status.init_SWE,_swe);

    // Computed index by hand (prior to dividing by SWE)
    const double index = 36.5924;

    // Check index computation
    EXPECT_NEAR(status.index,index/_swe,diff);

    // Check that the major melt count increments
    EXPECT_EQ(status.major_melt_count,1);

    // check that maximum per day is properly set
    EXPECT_NEAR(status.max_major_per_melt,index / infDays, diff);
    
    // Computation of inf is conditional. Most likely only one ius true at a time.
    EXPECT_NEAR(model.get_snow_inf(),expected_inf(yesterday_melt),diff);
};

TEST_F(CrackTest, ConstantInfCheck)
{
const double start_melt = _snowmelt * steps_per_day;
    status.daily_melt_total = start_melt;

    double yesterday;
    set_newday(true);
    for (int i = 0; i < steps_per_day; i++)
    {

        Crack model = run_a_step(_snowmelt);
        // checking that infiltration from hour i is equal to hour i-1 (skipping the first hour).
        if (i > 0)
            EXPECT_EQ(model.get_snow_inf(),yesterday) << "Current hour: " << i;
        set_newday(false);
        yesterday = model.get_snow_inf();
    }
};

TEST_F(CrackTest, MultiDayInfiltration)
{
    const double start_melt = 6.0;
    const std::vector<double> melt = {start_melt, 2.0 * start_melt, 3.0 * start_melt};
    for (int j = 0; j < melt.size(); j++)
    {
        status.daily_melt_total = melt[j];
        set_newday(true);
        double inf; 
        for (int i = 0; i < steps_per_day; i++)
        {
            Crack model = run_a_step(_snowmelt);
            set_newday(false);
                    
            inf = model.get_snow_inf();
            
            //Test infiltration on every hour over several major days
            EXPECT_NEAR(inf,expected_inf(melt[j]),diff);
            //inf and snow_inf should be equal
            EXPECT_EQ(inf,model.get_inf());
        }
    }
};

TEST_F(CrackTest, MajorMeltCounterTest)
{
    const double start_melt = 6.0;
    std::vector<double> melt;
    const int num_major = 12;
    for (int j = 0; j < num_major; j ++)
        melt.push_back(start_melt * (j + 1));
    
    int num_skip = 7;
    melt[num_skip] = 1.0;
    for (int j = 0; j < melt.size(); j++)
    {
        status.daily_melt_total = melt[j];
        set_newday(true);
        for (int i = 0; i < steps_per_day; i++)
        {
            Crack model = run_a_step(_snowmelt);
            set_newday(false);
                    
        }
        // One day has a non-major, expect major_melt_count to not increment
        if (j == num_skip)
            EXPECT_EQ(status.major_melt_count,j);
        else if (j + 2 < infDays) //Check if incremented
            EXPECT_EQ(status.major_melt_count,j+1);
        
        // checks that infiltration stops after the limit is found
        if (j > infDays + 1)
            EXPECT_EQ(status.current_inf,0.0) << "Day: " << j+1;
    }
};

TEST_F(CrackTest, RestrictedTestNewDay)
{
    double day_melt = 48.0;
    set_newday(true);
    double storage = 100.0;
    Crack model = restricted_or_unlimited(day_melt,storage);
    
    // Check that restricted moves all melt to runoff.
    EXPECT_DOUBLE_EQ(model.get_inf(),0.0);
    EXPECT_DOUBLE_EQ(model.get_snow_inf(),0.0);
    EXPECT_DOUBLE_EQ(model.get_runoff(),day_melt);
    EXPECT_DOUBLE_EQ(model.get_melt_runoff(),day_melt);

};

TEST_F(CrackTest,RestrictedTestManyDays)
{
    const double start_melt = 6.0;
    _soil_storage_at_freeze = 100.0;
    const std::vector<double> melt = {start_melt, 2.0 * start_melt, 3.0 * start_melt};
    for (int j = 0; j < melt.size(); j++)
    {
        status.daily_melt_total = melt[j];
        set_newday(true);
        double inf; 
        for (int i = 0; i < steps_per_day; i++)
        {
            Crack model = run_a_step(_snowmelt);
            set_newday(false);
                    
            inf = model.get_snow_inf();
            // Same as last, buyt many days.
            EXPECT_EQ(inf,0.0);
            EXPECT_EQ(model.get_snow_inf(),0.0);
            EXPECT_EQ(model.get_runoff(),melt[j]);
            EXPECT_EQ(model.get_melt_runoff(),melt[j]);

        }
    }
        
};

TEST_F(CrackTest,UnlimitedTestManyDays)
{
    const double start_melt = 6.0;
    _soil_storage_at_freeze = 0.0;
    const std::vector<double> melt = {start_melt, 2.0 * start_melt, 3.0 * start_melt};
    for (int j = 0; j < melt.size(); j++)
    {
        status.daily_melt_total = melt[j];
        set_newday(true);
        double inf; 
        for (int i = 0; i < steps_per_day; i++)
        {
            Crack model = run_a_step(_snowmelt);
            set_newday(false);
                    
            inf = model.get_inf();
            //Same as last but now everything infiltrates because it is unlimited conditions
            EXPECT_EQ(inf,melt[j]);
            EXPECT_EQ(model.get_snow_inf(),melt[j]);
            EXPECT_EQ(model.get_runoff(),0.0);
            EXPECT_EQ(model.get_melt_runoff(),0.0);

        }
    }
        
};

TEST_F(CrackTest,IceLensTest)
{
    set_newday(true);
    status.daily_melt_total = _snowmelt * steps_per_day;
    std::vector<double> T;

    for (int i = 0; i < steps_per_day; i++)
        T.push_back(-20.0);

    for (int i = 0; i < steps_per_day;i++)
    {
        _airtemp = T[i];
        Crack model = run_a_step(_snowmelt);
        set_newday(false);
        if (i > 0)
            EXPECT_EQ(status.tmax,T[0]); 
    }
    set_newday(true);
    Crack model = run_a_step(_snowmelt*3.0);

    // Check that the ice-lens is set by increasing major_melt_count safely above InfDays.
    EXPECT_EQ(status.major_melt_count,infDays+4);
    
};

TEST_F(CrackTest,RainOnSnowTest)
{
    status.init();
    status.begin_freeze();
    double rain_on_snow,sum,oldsum;
    sum = 0;
    
    for (int ii = 0; ii < 48; ++ii)
    {
        _rainfall = 2.0*ii;
        sum += _rainfall;
        _newday = ii % 24 == 0;
        if (_newday) 
        {   
            oldsum = sum;
            sum = 0.0;
        }
        Crack model = run_a_step(0.0);

        rain_on_snow = model.get_rain_on_snow();
        
        EXPECT_EQ(sum,status.daily_rain_total) << "Step: " << ii;
        if (_newday)
            EXPECT_EQ(oldsum,rain_on_snow) << "Step: " << ii;
    }
};

class CrackImplTest : public testing::Test
{
protected:

    CrackImplTest()
    {
        status.init();
    };
    CSVReader reader;
    Crack::info status;
    static constexpr double seconds_per_hour = 3600.0;
    double steps_per_day = 24;
    double rainfall;
    double snowfall;

    double major = 5.0;
    double min_swe_to_freeze = 25.0;
    unsigned int infDays = 6;
    bool AllowPriorInf = true;
    double lenstemp = -10.0;

    struct CRHM
    {
        double infil;
        double snowinfil;
        double melt_runoff;
        double runoff;
        double rain_on_snow;

        CRHM(const int& i,CSVReader& reader)
        {
            infil = reader.getValue<double>("infil",i);
            snowinfil = reader.getValue<double>("snowinfil",i) / 24;
            melt_runoff = reader.getValue<double>("meltrunoff",i) / 24;
            runoff = reader.getValue<double>("runoff",i);
            rain_on_snow = reader.getValue<double>("RainOnSnow",i);
        };
    };

    template<typename T>
    void print(std::string text,T val)
    {
        std::cout << text << val << std::endl;
    };
     

};

TEST_F(CrackImplTest,FullImplementTest)
{ 
    status.init();
    int start = 0;
    int end = 140000;
    for (int i = start; i < end; ++i)
    { 
        //std::cout << " " <<std::endl;
        print("TIME STEP: ", i); 
        double rainfall = reader.getValue<double>("net_rain",i);
        double snowmelt = reader.getValue<double>("snowmeltD",i) / 24;
        double swe = reader.getValue<double>("SWE",i);
        double soil_storage_at_freeze = 50;
        double airtemp = reader.getValue<double>("hru_t",i);
        bool crackon = reader.getValue<bool>("crackon",i);
        std::string datetime = reader.getValue<std::string>("datetime",i);
        //std::cout << datetime << std::endl;
        CRHM crhm(i,reader);
        print("rainfall: ", rainfall);
        print("snowmelt: ", snowmelt);
        print("swe: ", swe);

        double runoff = 0.0;
        double melt_runoff = 0.0;
        double inf = 0.0;
        double snowinf = 0.0;
        double rain_on_snow = 0.0;
        bool is_day_over = i % 24 == 23;
        print("Day over tracker: ", i % 24);

        if ((swe > min_swe_to_freeze && !status.frozen && is_day_over) || (crackon && i == start) )
        {
            status.begin_freeze();
            status.end_freeze_tomorrow = false;
        }
        
        print("snowinf: ", snowinf);
        print("inf: ", inf);
        print("melt_runoff: ", melt_runoff);
        print("runoff: ", runoff);
        print("frozen: ",status.frozen);
        print("major melt count: ",status.major_melt_count);
        print("Total input: ", snowmelt+rainfall);
        
        if (status.frozen)
        {
            Crack crack(major, min_swe_to_freeze, infDays, 
                    AllowPriorInf, lenstemp,steps_per_day,status);
            
            crack.init_inputs(snowmelt, rainfall, swe, soil_storage_at_freeze,
                    airtemp, is_day_over); 
            status.daily_melt_total = snowmelt * steps_per_day;
            crack.is_CRHM_compare_test = true;
            crack.run();

            runoff = crack.get_runoff() / steps_per_day;
            melt_runoff = crack.get_melt_runoff() / steps_per_day;
            inf = crack.get_inf() / steps_per_day;
            snowinf = crack.get_snow_inf() / steps_per_day;
            rain_on_snow = crack.get_rain_on_snow();
        
            if (is_day_over && swe <= 0.0 && status.major_melt_count > 0)
                status.end_freeze();
            //if (swe <= 0.0 && status.major_melt_count > 0)
            //    status.end_freeze_tomorrow = true;

        }
        else
        {
            crhm.infil = 0.0;
            crhm.snowinfil = 0.0;
            crhm.melt_runoff = 0.0;
            crhm.runoff = 0.0;
            crhm.rain_on_snow = 0.0;
        };
        print("Melt total: ",status.daily_melt_total);
        print("Rain total: ",status.daily_rain_total);
        print("index: ", status.index);
        print("Max major per melt: ", status.max_major_per_melt);
        print("init_SWE", status.init_SWE);

        print("snowinf: ", snowinf);
        print("inf: ", inf);
        print("melt_runoff: ", melt_runoff);
        print("runoff: ", runoff);
        double diff = 1e-5;
        EXPECT_NEAR(crhm.infil,inf - snowinf,diff) << "Step: " << i;
        EXPECT_NEAR(crhm.snowinfil,snowinf,diff) << "Step: " << i;
        EXPECT_NEAR(crhm.melt_runoff,melt_runoff,diff) << "Step: " << i;
        EXPECT_NEAR(crhm.runoff,runoff - melt_runoff,diff) << "Step: " << i;
        EXPECT_EQ(status.frozen,crackon) << "Step: " << i;
        EXPECT_NEAR(crhm.rain_on_snow,rain_on_snow,diff) << "Step: " << i;

        
    }
};
