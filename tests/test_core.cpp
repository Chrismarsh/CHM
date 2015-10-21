
#include "../src/core.h"
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

TEST_F(CoreTest,ThrowsOnInvalidFile)
{


    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "invalid_path.config"
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
                    (char*) "test_config_file.config"
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
                    (char*) "test_missing_modules.config"
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
                    (char*) "test_empty_modules.config"
            };

    ASSERT_ANY_THROW(c1.init(3,argv));
}

//Tests proper handling when an optional section is missing
TEST_F(CoreTest,OptionalSectionDebug)
{
    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "test_optional_section.config"
            };

    ASSERT_NO_THROW(c1.init(3,argv));
}

//Tests proper handling when an optional section is missing
TEST_F(CoreTest,BrokenConfig)
{
    core c1;
    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-f",
                    (char*) "test_malformed_config.config"
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
                    (char*) "test_config_file.config",
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
                    (char*) "test_config_file.config",
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
                    (char*) "test_config_file.config",
                    (char*) "-c",
                    (char*) "config.Harder_precip_phase.const.b:3.14",
                    (char*) "-c",
                    (char*) "nproc:2"

            };


    ASSERT_NO_THROW(c1.init(6,argv));

    ASSERT_EQ(c1._cfg.get<double>("config.Harder_precip_phase.const.b"),3.14);
    ASSERT_EQ(c1._cfg.get<int>("nproc"),2);

}

//tests positional argument to pass config without -f, but with the arguments in a funny order
TEST_F(CoreTest,CmdlOptionsPositional2)
{
    core c1;

    char* argv[] =
            {
                    (char*) "/path/to/CHM",
                    (char*) "-c",
                    (char*) "config.Harder_precip_phase.const.b:3.14",
                    (char*) "test_config_file.config",
                    (char*) "-c",
                    (char*) "nproc:2"
            };


    ASSERT_NO_THROW(c1.init(6,argv));

    ASSERT_EQ(c1._cfg.get<double>("config.Harder_precip_phase.const.b"),3.14);
    ASSERT_EQ(c1._cfg.get<int>("nproc"),2);

}
