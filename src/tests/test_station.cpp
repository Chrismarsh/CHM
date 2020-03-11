//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//


#include "station.hpp"
#include "gtest/gtest.h"

class StationTest : public testing::Test
{
protected:


    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);

    }

    station s0;
    station s1;
    std::set<std::string> vars={"t","rh","u"};
};


//tests for proper default init
TEST_F(StationTest, DefaultInit)
{
    EXPECT_EQ(0, s0.x());
    EXPECT_EQ(0, s0.y());
    EXPECT_EQ(0, s0.z());
}


TEST_F(StationTest, Equal)
{
    station s1("upper clearing", 60.52163, -135.197151, -1305,vars );
    station s2("granger", 60.52163, -135.197151, -1305, vars );

    station s3("upper clearing", 60.52163, -135.197151, -1305, vars);


    EXPECT_FALSE(s1==s2);
    EXPECT_TRUE(s1==s3);


}