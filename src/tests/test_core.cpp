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


#include "core.hpp"
#include "gtest/gtest.h"
#include <stdlib.h>
#include <string>
#include <utility>
/**
 * Tests various functions of the core object including:
 * - Reading config files
 *
 */


class CoreTest : public testing::Test
{

protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);

//        ASSERT_NO_THROW(c0.init(argc,argv));
    }
   // core c0;

};

//TEST_F(CoreTest,CmdlHelp)
//{
//
//
//    core c1;
//    char* argv[] =
//            {
//                    (char*) "/path/to/CHM",
//                    (char*) "-h"
//            };
//    ASSERT_ANY_THROW(c1.init(2,argv));
//}
//
//TEST_F(CoreTest,CmdlVersion)
//{
//
//
//    core c1;
//    char* argv[] =
//            {
//                    (char*) "/path/to/CHM",
//                    (char*) "-v"
//            };
//    ASSERT_ANY_THROW(c1.init(2,argv));
//}

TEST_F(CoreTest,ThrowsOnInvalidFile)
{


    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "invalid_path.json"
            };
    ASSERT_ANY_THROW(c1.init(3,argv));
}

TEST_F(CoreTest,NoThrowOnValidFile)
{
    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "test_config_file.json"
            };

    ASSERT_NO_THROW(c1.init(3,argv));
}

//Tests proper handling when missing a required section: Modules
TEST_F(CoreTest,MissingSectionModules)
{
    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "test_config_file_missing_modules.json"
            };


    ASSERT_ANY_THROW(c1.init(3,argv));
}

//Tests proper handling when a required section is empty: Modules
TEST_F(CoreTest,EmptySectionModules)
{
    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "test_config_file_empty_modules.json"
            };

    ASSERT_ANY_THROW(c1.init(3,argv));
}

////Tests proper handling when an optional section is missing
//TEST_F(CoreTest,OptionalSectionDebug)
//{
//    core c1;
//    char* argv[] =
//            {
//                    (char*) "/path/to/CHM",
//                    (char*) "-f",
//                    (char*) "test_optional_section.json"
//            };
//
//    ASSERT_NO_THROW(c1.init(3,argv));
//}

//Tests proper handling when an optional section is missing
TEST_F(CoreTest,BrokenConfig)
{
    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "test_config_file_malformed.json"
            };

    ASSERT_ANY_THROW(c1.init(3,argv));
}

TEST_F(CoreTest,CmdlOptions)
{
    core c1;

    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "test_config_file.json",
                    (char*) "-c",
                    (char*) "config.Harder_precip_phase.const.b:3.14",
                    (char*) "-c",
                    (char*) "nproc:2"

            };

    ASSERT_NO_THROW(c1.init(7,argv));

    ASSERT_EQ(c1._cfg.get<double>("config.Harder_precip_phase.const.b"),3.14);
    ASSERT_EQ(c1._cfg.get<int>("nproc"),2);

}

//same as above but tests -c being used twice
TEST_F(CoreTest,CmdlOptions2)
{
    core c1;

    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "test_config_file.json",
                    (char*) "-c",
                    (char*) "config.Harder_precip_phase.const.b:3.14",
                    (char*) "-c",
                    (char*) "nproc:2"

            };


    ASSERT_NO_THROW(c1.init(6,argv));

    ASSERT_EQ(c1._cfg.get<double>("config.Harder_precip_phase.const.b"),3.14);
    ASSERT_EQ(c1._cfg.get<int>("nproc"),2);

}

//tests positional argument to pass config without -f
TEST_F(CoreTest,CmdlOptionsPositional)
{
    core c1;

    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "test_config_file.json",
            };


    ASSERT_NO_THROW(c1.init(2,argv));


}

//tests positional argument to pass config without -f, but with the arguments in a funny order
TEST_F(CoreTest,CmdlOptionsPositional2)
{
    core c1;

    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "test_config_file.json",
                    (char*) "-c",
                    (char*) "config.Harder_precip_phase.const.b:3.14",
                    (char*) "-c",
                    (char*) "nproc:2"
            };


    ASSERT_NO_THROW(c1.init(6,argv));

    ASSERT_EQ(c1._cfg.get<double>("config.Harder_precip_phase.const.b"),3.14);
    ASSERT_EQ(c1._cfg.get<int>("nproc"),2);

}

//tests erasing keys
TEST_F(CoreTest,CmdlOptionsErase)
{
    core c1;

    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-c",
                    (char*) "config.Harder_precip_phase.const.b:3.14",
                    (char*) "-f",
                    (char*) "test_config_file.json",
                    (char*) "-c",
                    (char*) "nproc:2",
                    (char*) "-r",
                    (char*) "nproc",
                    (char*) "-r",
                    (char*) "output",
                    (char*) "-r",
                    (char*) "global.UTC_offset"

            };



    ASSERT_NO_THROW(c1.init(13,argv));
    ASSERT_EQ(c1._cfg.get<double>("config.Harder_precip_phase.const.b"),3.14);

    //removals are done after -c so, this should not exist

    ASSERT_ANY_THROW(c1._cfg.get<int>("nproc"));
    ASSERT_ANY_THROW(c1._cfg.get<int>("global.UTC_offset"));
    ASSERT_ANY_THROW(c1._cfg.get_child("output"));  //confirms the entire section got nuked


}
