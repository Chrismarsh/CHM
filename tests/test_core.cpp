
#include "../src/core.h"
#include "gtest/gtest.h"

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
            
        ASSERT_NO_THROW(c0.read_config_file("test_config_file.config"));
    }
    core c0;

};

TEST_F(CoreTest,ThrowsOnInvalidFile)
{
    core c1;
    ASSERT_ANY_THROW(c1.read_config_file("invalid_path.config"));
}

TEST_F(CoreTest,NoThrowOnValidFile)
{
    core c1;
    ASSERT_NO_THROW(c1.read_config_file("test_config_file.config"));
}

//Tests proper handling when missing a required section: Modules
TEST_F(CoreTest,MissingSectionModules)
{
    core c1;
    ASSERT_ANY_THROW(c1.read_config_file("test_missing_modules.config"));
}

//Tests proper handling when a required section is empty: Modules
TEST_F(CoreTest,EmptySectionModules)
{
    core c1;
    ASSERT_ANY_THROW(c1.read_config_file("test_empty_modules.config"));
}

//Tests proper handling when an optional section is missing
TEST_F(CoreTest,OptionalSectionDebug)
{
    core c1;
    ASSERT_NO_THROW(c1.read_config_file("test_optional_section.config"));
}

//Tests proper handling when an optional section is missing
TEST_F(CoreTest,BrokenConfig)
{
    core c1;
    ASSERT_ANY_THROW(c1.read_config_file("test_malformed_config.config"));
}