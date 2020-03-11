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

#include "metdata.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <string>
#include <algorithm>

class MetdataTest : public testing::Test
{
  protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);

    }
    std::string proj4str = "+proj=utm +zone=8 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs ";

};

//basic default init sanity checks
TEST_F(MetdataTest, LoadAsciiData)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude =  -135.184652;
    station.elevation = 1559;

    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s,-8));

}

TEST_F(MetdataTest, ASCII_TwoStationDuplicatedID)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude =  -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station2;
    station2.path = "test_met_data_longer2.txt";
    station2.latitude = 60.56726;
    station2.longitude =  -135.184652;
    station2.elevation = 1559;
    station2.id = "station1";

    std::vector<metdata::ascii_metdata> s;

    s.push_back(station2);
    s.push_back(station);

    ASSERT_ANY_THROW(md.load_from_ascii(s,-8));

}

TEST_F(MetdataTest, ASCII_TwoStation)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude =  -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station2;
    station2.path = "test_met_data_longer2.txt";
    station2.latitude = 60.56726;
    station2.longitude =  -135.184652;
    station2.elevation = 1559;
    station2.id = "station2";

    std::vector<metdata::ascii_metdata> s;

    s.push_back(station2);
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s,-8));

}
TEST_F(MetdataTest, ASCII_StartEndTime)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude =  -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);


    ASSERT_NO_THROW(md.load_from_ascii(s,-8));

    auto start_time = md.start_time();
    auto end_time = md.end_time();

    auto s_true = boost::posix_time::from_iso_string("20101001T090000");
    auto e_true = boost::posix_time::from_iso_string("20101001T120000");
    ASSERT_EQ(start_time,s_true);
    ASSERT_EQ(end_time,e_true);

}


TEST_F(MetdataTest, ASCII_TwoStationStartEndTime)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude =  -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station2;
    station2.path = "test_met_data_longer2.txt";
    station2.latitude = 60.56726;
    station2.longitude =  -135.184652;
    station2.elevation = 1559;
    station2.id = "station2";

    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);
    s.push_back(station2);

    ASSERT_NO_THROW(md.load_from_ascii(s,-8));

    boost::posix_time::ptime start_time,end_time;
    std::tie(start_time, end_time) = md.start_end_time();

    auto s_true = boost::posix_time::from_iso_string("20101001T100000");
    auto e = boost::posix_time::from_iso_string("20101001T120000");
    ASSERT_EQ(start_time,s_true);
    ASSERT_EQ(end_time,e);

}
TEST_F(MetdataTest, ASCII_TestCurrentTimeStr)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));

    auto cur_time = md.current_time_str();

    ASSERT_EQ(cur_time,"20101001T090000");
}

TEST_F(MetdataTest, ASCII_TestNext)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station2;
    station2.path = "test_met_data_longer2.txt";
    station2.latitude = 60.56726;
    station2.longitude = -135.184652;
    station2.elevation = 1559;
    station2.id = "station2";

    std::vector<metdata::ascii_metdata> s;

    s.push_back(station2);
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));
    md.next(); // load first timestep

    std::string cur_time = "";
    //returns false when no ts left


    // First Ts = 20101001T100000
    cur_time = md.current_time_str();
    ASSERT_EQ(cur_time,"20101001T100000");

    // 20101001T110000
    bool done = !md.next();
    ASSERT_EQ(done,false);
    cur_time = md.current_time_str();
    ASSERT_EQ(cur_time,"20101001T110000");

    // 20101001T120000
    done = !md.next();
    ASSERT_EQ(done,false);
    cur_time = md.current_time_str();
    ASSERT_EQ(cur_time,"20101001T120000");
}

TEST_F(MetdataTest, ASCII_TestAccessData)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));
    md.next(); // load first timestep

    ASSERT_DOUBLE_EQ(md.at(0)->operator[]("Qsi"_s),112.101);

}

TEST_F(MetdataTest, ASCII_TestStartTimeStr)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));

    ASSERT_EQ( md.start_time_str(),"20101001T090000");

}

TEST_F(MetdataTest, ASCII_TestEndTimeStr)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));

    ASSERT_EQ( md.end_time_str(),"20101001T120000");

}

TEST_F(MetdataTest, ASCII_TestNStations)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station2;
    station2.path = "test_met_data_longer2.txt";
    station2.latitude = 60.56726;
    station2.longitude = -135.184652;
    station2.elevation = 1559;
    station2.id = "station2";

    std::vector<metdata::ascii_metdata> s;

    s.push_back(station2);
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));

    ASSERT_EQ( md.nstations(),2);

}

TEST_F(MetdataTest, ASCII_TestDt)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station2;
    station2.path = "test_met_data_longer2.txt";
    station2.latitude = 60.56726;
    station2.longitude = -135.184652;
    station2.elevation = 1559;
    station2.id = "station2";

    std::vector<metdata::ascii_metdata> s;

    s.push_back(station2);
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));

    ASSERT_EQ( md.dt_seconds(),3600);

}

TEST_F(MetdataTest, ASCII_TestMissingTs)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "missing_timestep.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";


    std::vector<metdata::ascii_metdata> s;

    s.push_back(station);

    ASSERT_ANY_THROW(md.load_from_ascii(s, -8));

}

TEST_F(MetdataTest, ASCII_TestInconsistentTs)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station3;
    station3.path = "test_met_data_longer3_3hr_dt.txt";
    station3.latitude = 60.56726;
    station3.longitude = -135.184652;
    station3.elevation = 1559;
    station3.id = "station3";


    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);
    s.push_back(station3);

    ASSERT_ANY_THROW(md.load_from_ascii(s, -8));

}

TEST_F(MetdataTest, ASCII_TestSubset)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";

    metdata::ascii_metdata station2;
    station2.path = "test_met_data_longer2.txt";
    station2.latitude = 60.56726;
    station2.longitude = -135.184652;
    station2.elevation = 1559;
    station2.id = "station2";


    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);
    s.push_back(station2);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));

    boost::posix_time::ptime start = boost::posix_time::from_iso_string("20101001T110000");
    boost::posix_time::ptime end = boost::posix_time::from_iso_string("20101001T120000");

    ASSERT_NO_THROW(md.subset(start,end));

    ASSERT_EQ(md.start_time_str(),"20101001T110000");
    ASSERT_EQ(md.end_time_str(),"20101001T120000");
    ASSERT_EQ(md.n_timestep(),2);

}

TEST_F(MetdataTest, ASCII_TestListVars)
{
    metdata md(proj4str);

    metdata::ascii_metdata station;
    station.path = "test_met_data_longer1.txt";
    station.latitude = 60.56726;
    station.longitude = -135.184652;
    station.elevation = 1559;
    station.id = "station1";


    std::vector<metdata::ascii_metdata> s;
    s.push_back(station);

    ASSERT_NO_THROW(md.load_from_ascii(s, -8));

    auto vars = md.list_variables();
    std::set<std::string> variable_list = {"t",   "rh",  "u",   "p",   "Qsi","T_g"};

    ASSERT_EQ(vars.size(),variable_list.size());
    ASSERT_TRUE(vars == variable_list );

}



TEST_F(MetdataTest, NC_TestBasicLoad)
{
    metdata md(proj4str);

    ASSERT_NO_THROW(md.load_from_netcdf("GEM-CHM_2p5_snowcast_2018011506_2018011605.nc"));

    md.next(); // load first timestep

    ASSERT_EQ(md.nstations(), 151*151 );

    auto vars = md.list_variables();
    std::set<std::string> variable_list = {"p","press","Qli","Qsi","Qsi_diff","rh","t","u","vw_dir"};

    ASSERT_EQ(vars.size(),variable_list.size());
    ASSERT_TRUE(vars == variable_list );

    ASSERT_EQ(md.start_time_str(),"20180115T060000");
    ASSERT_EQ(md.end_time_str(),"20180116T050000");

    ASSERT_EQ(md.current_time_str(),"20180115T060000");
}

TEST_F(MetdataTest, NC_TestAcess)
{
    metdata md(proj4str);

    ASSERT_NO_THROW(md.load_from_netcdf("GEM-CHM_2p5_snowcast_2018011506_2018011605.nc"));

    md.next(); // load first timestep

    auto value = (*md.at(0))["t"];
    ASSERT_DOUBLE_EQ(value, -1.4661178588867188);

    value = (*md.at(150+150*151))["t"];
    ASSERT_DOUBLE_EQ(value, -18.4973678588867188);

}

TEST_F(MetdataTest, NC_TestNext)
{
    metdata md(proj4str);

    ASSERT_NO_THROW(md.load_from_netcdf("GEM-CHM_2p5_snowcast_2018011506_2018011605.nc"));
    md.next(); // load first timestep


    auto value = (*md.at(0))["t"];
    ASSERT_DOUBLE_EQ(value, -1.4661178588867188);

    ASSERT_NO_THROW(md.next());

    ASSERT_EQ(md.current_time_str(),"20180115T070000");

    value = (*md.at(0))["t"];
    ASSERT_DOUBLE_EQ(value, -1.5229721069335938);

    value = (*md.at(150+150*151))["t"];
    ASSERT_DOUBLE_EQ(value, -18.5444564819335938);

}


TEST_F(MetdataTest, NC_TestnTimeSteps)
{
    metdata md(proj4str);

    ASSERT_NO_THROW(md.load_from_netcdf("GEM-CHM_2p5_snowcast_2018011506_2018011605.nc"));

    ASSERT_EQ(md.n_timestep(),24);
}

TEST_F(MetdataTest, NC_TestPrune)
{
    metdata md(proj4str);

    ASSERT_NO_THROW(md.load_from_netcdf("GEM-CHM_2p5_snowcast_2018011506_2018011605.nc"));

    std::unordered_set<std::string> ids_to_remove;

    for(size_t i = 1; i < 151; i++)
    {
        for(size_t j = 0; j < 151; j++)
        {
            size_t index = i + j * 151;
            std::string station_name = std::to_string(index); // these don't really have name
            ids_to_remove.insert(station_name);
        }
    }

    md.prune_stations(ids_to_remove);

    ASSERT_EQ((151*151)-(151*150),md.nstations());
    ASSERT_EQ(md.stations().at(0)->ID(),"0");
}