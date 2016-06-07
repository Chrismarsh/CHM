
#include "../timeseries/daily.hpp"
#include "../timeseries/timeseries.hpp"
#include "gtest/gtest.h"

class DailyTest : public testing::Test
{
protected:
    timeseries ts;

    virtual void SetUp()
    {

        logging::core::get()->set_logging_enabled(false);

        ASSERT_NO_THROW(ts.open("test_daily_vv_oct2010.txt"));
    }


    ~DailyTest()
    {

    }
};

TEST_F(DailyTest,StartofDay)
{

    auto itr = ts.begin();
    itr ++;

    auto sod = daily::start_of_day(ts,itr);

    ASSERT_EQ(2,sod->day());
    ASSERT_EQ(0,sod->hour());
    ASSERT_EQ(0,sod->min());
    ASSERT_DOUBLE_EQ(9.67,sod->get("t"));

    itr = ts.begin();
    itr=itr+32;
    ASSERT_DOUBLE_EQ(6.524, itr->get("t")); //absolutely make sure we're where we think we are. Which is 20101003T080000

    sod = daily::start_of_day(ts,itr);
    ASSERT_DOUBLE_EQ(8.635, sod->get("t"));
    ASSERT_EQ(3,sod->day());
    ASSERT_EQ(0,sod->hour());
    ASSERT_EQ(0,sod->min());
}

TEST_F(DailyTest,EndofDay)
{

    auto itr = ts.begin();
    itr ++;

    auto sod = daily::end_of_day(ts,itr);

    ASSERT_EQ(2,sod->day());
    ASSERT_EQ(23,sod->hour());
    ASSERT_EQ(0,sod->min());
    ASSERT_DOUBLE_EQ(9.990,sod->get("t"));

    itr = ts.begin();
    itr=itr+32;
    ASSERT_DOUBLE_EQ(6.524, itr->get("t")); //absolutely make sure we're where we think we are. Which is 20101003T080000


    sod = daily::end_of_day(ts,itr);
    ASSERT_DOUBLE_EQ(10.550, sod->get("t"));
    ASSERT_EQ(3,sod->day());
    ASSERT_EQ(23,sod->hour());
    ASSERT_EQ(0,sod->min());
}

TEST_F(DailyTest,DailyMax)
{

    double max[30] = {20.05,17.265,14.828,12.885,18.04,17.743,13.053,12.47,14.785,5.637,5.223,11.05,13.231,8.486,0.464,3.817,6.783,8.249,10.453,12.419,10.243,7.011,4.178,1.352,-0.04,3.599,0.404,2.574,4.478,4.859};

    for(auto itr=ts.begin(); itr != ts.end(); itr++)
    {
        double m = daily::max(ts, itr, "t");
        ASSERT_DOUBLE_EQ(m,max[itr->day()-2]);
    }


}

TEST_F(DailyTest,DailyMin)
{

    double min[30] = {7.612, 6.524, 6.85, 0.874, 3.462, 5.391, 5.378, 5.014, 5.242, 0.05, -0.46, 0.898, 4.162, -4.272, -6.116, -5.756, -2.273, 1.43, 1.345, 2.532, 0.408, 0.379, 0.349, -1.769, -2.484, -4.185, -5.923, -3.27, -4.337, -2.789};

    for (auto itr = ts.begin(); itr != ts.end(); itr++)
    {
        double m = daily::min(ts, itr, "t");

        ASSERT_DOUBLE_EQ(m, min[itr->day() - 2]);

    }
}


