

#include "timeseries.hpp"
#include "gtest/gtest.h"


class TimeseriesTest : public testing::Test
{
protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
        ASSERT_NO_THROW(s1.open("test_met_data.txt"));
    }

    timeseries s0;
    timeseries s1;
};


//tests for proper default init

TEST_F(TimeseriesTest, DefaultInit)
{
    EXPECT_FALSE(s0.is_open());
}

//Test that a file that exist doesn't throw
//Test that the internal iterator doesn't throw on dereference

TEST_F(TimeseriesTest, ExistFileOpens)
{
    ASSERT_NO_THROW(s0.open("test_met_data.txt"));
}

//Test that opening a file that doesn't exist throws

TEST_F(TimeseriesTest, NoExistFileOpens)
{
    ASSERT_ANY_THROW(s0.open("file_that_doesnt_exist.txt"));
}
TEST_F(TimeseriesTest, init_new_variable)
{
    timeseries s;
    s.open("test_met_data.txt");

    s.init_new_variable("new_var");

    ASSERT_DOUBLE_EQ(-9999,s.at("new_var",0));
    s.at("new_var",0) = 42.314159;
    ASSERT_DOUBLE_EQ(42.314159,s.at("new_var",0));

    auto itr = s1.begin();
    ASSERT_DOUBLE_EQ(-8.1, itr->get("t"));
}

TEST_F(TimeseriesTest, read_comma_delimited)
{
    timeseries s;
    s.open("vv_oct12010_ts1_comma_delim.txt");

    auto itr = s.begin();

    ASSERT_DOUBLE_EQ(3.77, itr->get("t"));
    ASSERT_DOUBLE_EQ(70.720, itr->get("rh"));
    ASSERT_DOUBLE_EQ(1.444, itr->get("u"));
    ASSERT_DOUBLE_EQ(0.0000, itr->get("p"));
    ASSERT_DOUBLE_EQ(744.954, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(6.326, itr->get("T_g"));

    ASSERT_EQ(02, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());

    ASSERT_EQ(13, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());
}

TEST_F(TimeseriesTest, Advance)
{
    auto itr = s1.begin();
    itr++;

    ASSERT_DOUBLE_EQ(14.268, itr->get("t"));
    ASSERT_DOUBLE_EQ(21.739, itr->get("rh"));
    ASSERT_DOUBLE_EQ(1.850, itr->get("u"));
    ASSERT_DOUBLE_EQ(0.0000, itr->get("p"));
    ASSERT_DOUBLE_EQ(650.219, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(4.687, itr->get("T_g"));

    ASSERT_EQ(01, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());

    ASSERT_EQ(13, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());
}
//ensure we read in the right thing

TEST_F(TimeseriesTest, FileContents)
{
    auto itr = s1.begin();
    ASSERT_DOUBLE_EQ(-8.1, itr->get("t"));
    ASSERT_DOUBLE_EQ(40.282, itr->get("rh"));
    ASSERT_DOUBLE_EQ(0.647, itr->get("u"));
    ASSERT_DOUBLE_EQ(0.1031, itr->get("p"));
    ASSERT_DOUBLE_EQ(112.101, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(3.859, itr->get("T_g"));

    ASSERT_EQ(01, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());
    
    ASSERT_EQ(9, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());

    itr++;

    ASSERT_DOUBLE_EQ(14.268, itr->get("t"));
    ASSERT_DOUBLE_EQ(21.739, itr->get("rh"));
    ASSERT_DOUBLE_EQ(1.850, itr->get("u"));
    ASSERT_DOUBLE_EQ(0.0000, itr->get("p"));
    ASSERT_DOUBLE_EQ(650.219, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(4.687, itr->get("T_g"));

    ASSERT_EQ(01, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());
    
    ASSERT_EQ(13, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());

}

//Tests if we can change the timeseries
TEST_F(TimeseriesTest, SetContents)
{
    auto itr = s1.begin();
    
    itr->set("t",11.43);
    itr->set("p",3.0);
    itr++;
    itr->set("Qsi",743.32);
    itr->set("T_g",2.4);
    
    itr = s1.begin();
    ASSERT_DOUBLE_EQ(11.43, itr->get("t"));
    ASSERT_DOUBLE_EQ(40.282, itr->get("rh"));
    ASSERT_DOUBLE_EQ(0.647, itr->get("u"));
    ASSERT_DOUBLE_EQ(3.0, itr->get("p"));
    ASSERT_DOUBLE_EQ(112.101, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(3.859, itr->get("T_g"));

    ASSERT_EQ(01, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());
    
    ASSERT_EQ(9, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());

    itr++;

    ASSERT_DOUBLE_EQ(14.268, itr->get("t"));
    ASSERT_DOUBLE_EQ(21.739, itr->get("rh"));
    ASSERT_DOUBLE_EQ(1.850, itr->get("u"));
    ASSERT_DOUBLE_EQ(0.0000, itr->get("p"));
    ASSERT_DOUBLE_EQ(743.32, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(2.4, itr->get("T_g"));

    ASSERT_EQ(01, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());
    
    ASSERT_EQ(13, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());

}


//tests access of the underlying iterators
TEST_F(TimeseriesTest, GetItr)
{
    auto itr = s1.begin()->get_itr("t");
    
    ASSERT_DOUBLE_EQ(*itr,-8.100);
    itr++;
    ASSERT_DOUBLE_EQ(*itr,14.268);
    
    itr = s1.begin()->get_itr("rh");
    ASSERT_DOUBLE_EQ(*itr,40.282);
    itr++;
    ASSERT_DOUBLE_EQ(*itr,21.739);
    
    itr = s1.begin()->get_itr("u");
    ASSERT_DOUBLE_EQ(*itr,0.647);
    itr++;
    ASSERT_DOUBLE_EQ(*itr,1.85);
    
    itr = s1.begin()->get_itr("p");
    ASSERT_DOUBLE_EQ(*itr,0.1031);
    itr++;
    ASSERT_DOUBLE_EQ(*itr,0.0);
    
    itr = s1.begin()->get_itr("Qsi");
    ASSERT_DOUBLE_EQ(*itr,112.101);
    itr++;
    ASSERT_DOUBLE_EQ(*itr,650.219);
    
    itr = s1.begin()->get_itr("T_g");
    ASSERT_DOUBLE_EQ(*itr,3.859);
    itr++;
    ASSERT_DOUBLE_EQ(*itr,4.687);

}

//test we are returning the correct length of the timeseries
TEST_F(TimeseriesTest, Length)
{
    EXPECT_EQ(2, s1.get_timeseries_length());
}

TEST_F(TimeseriesTest, ToFile)
{
    s1.to_file("test_s1_output.txt");
    timeseries input;
    ASSERT_NO_THROW(input.open("test_s1_output.txt"));
    
    auto itr = s1.begin();
    ASSERT_DOUBLE_EQ(-8.1, itr->get("t"));
    ASSERT_DOUBLE_EQ(40.282, itr->get("rh"));
    ASSERT_DOUBLE_EQ(0.647, itr->get("u"));
    ASSERT_DOUBLE_EQ(0.1031, itr->get("p"));
    ASSERT_DOUBLE_EQ(112.101, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(3.859, itr->get("T_g"));

    ASSERT_EQ(01, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());
    
    ASSERT_EQ(9, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());

    itr++;

    ASSERT_DOUBLE_EQ(14.268, itr->get("t"));
    ASSERT_DOUBLE_EQ(21.739, itr->get("rh"));
    ASSERT_DOUBLE_EQ(1.850, itr->get("u"));
    ASSERT_DOUBLE_EQ(0.0000, itr->get("p"));
    ASSERT_DOUBLE_EQ(650.219, itr->get("Qsi"));
    ASSERT_DOUBLE_EQ(4.687, itr->get("T_g"));

    ASSERT_EQ(01, itr->day());
    ASSERT_EQ(10, itr->month());
    ASSERT_EQ(2010, itr->year());
    
    ASSERT_EQ(13, itr->hour());
    ASSERT_EQ(0, itr->min());
    ASSERT_EQ(0, itr->sec());
}

TEST_F(TimeseriesTest, MissingTimeStep)
{
    ASSERT_ANY_THROW(s0.open("missing_timestep.txt"));
}