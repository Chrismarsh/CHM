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

#include "netcdf.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <algorithm>

class NetCDFTest : public testing::Test
{
  protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);

        nc.open_GEM("GEM-CHM_2p5_snowcast_2018011506_2018011605.nc");

    }

    netcdf nc;

};


TEST_F(NetCDFTest, VariableNames)
{
    auto vars = nc.get_variable_names();
    std::set<std::string> variable_list = {"p","press","Qli","Qsi","Qsi_diff","rh","t","u","vw_dir"};

    ASSERT_EQ(vars.size(),variable_list.size());
    ASSERT_TRUE(vars == variable_list );

}

TEST_F(NetCDFTest, CoordinateNames)
{
    auto vars = nc.get_coordinate_names();
    std::set<std::string> variable_list = {"datetime", "xgrid_0", "ygrid_0",};

    ASSERT_EQ(vars.size(),variable_list.size());

    ASSERT_TRUE(vars == variable_list );

}

TEST_F(NetCDFTest, StartTime)
{
    auto start = nc.get_start();

    ASSERT_EQ("20180115T060000",boost::posix_time::to_iso_string(start));

}

TEST_F(NetCDFTest, EndTime)
{
    auto end = nc.get_end();

    ASSERT_EQ("20180116T050000",boost::posix_time::to_iso_string(end));

}

TEST_F(NetCDFTest, Dt)
{
    auto dt = nc.get_dt();

    ASSERT_EQ(3600, dt.total_seconds());

}

TEST_F(NetCDFTest, ntimesteps)
{
    auto nts = nc.get_ntimesteps();

    ASSERT_EQ(24, nts);

}

TEST_F(NetCDFTest, access)
{
    auto time = boost::posix_time::from_iso_string("20180115T060000");
    auto value = nc.get_var("t",time,0,0);
    ASSERT_DOUBLE_EQ(value, -1.4661178588867188);

    value = nc.get_var("t",time,150,150);
    ASSERT_DOUBLE_EQ(value, -18.4973678588867188);

}

TEST_F(NetCDFTest, access_time_offset)
{
    auto time = boost::posix_time::from_iso_string("20180116T050000"); // last ts
    auto value = nc.get_var("t",time,0,0);
    ASSERT_DOUBLE_EQ(value, 0.3009796142578125);

    value = nc.get_var("t",time,150,150);
    ASSERT_DOUBLE_EQ(value, -11.3069305419921875);

}