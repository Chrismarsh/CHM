
#include "../src/station.hpp"
#include "gtest/gtest.h"

TEST(station, DefaultCtor)
{
    station s;
    EXPECT_EQ(0,s.get_x());
    EXPECT_EQ(0,s.get_y());
    EXPECT_EQ(0,s.get_z());
}