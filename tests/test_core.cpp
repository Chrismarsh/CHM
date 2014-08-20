
#include "../src/core.h"
#include "gtest/gtest.h"

class CoreTest : public testing::Test
{
protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
            
        ASSERT_NO_THROW(c0.read_config_file("test_config_file.config"));
    }
    core c0;
    void ReadConfigFile_Debug();
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

void CoreTest::ReadConfigFile_Debug()
{
    ASSERT_EQ(c0._log_level,debug);
}
//test different sections
TEST_F(CoreTest,ReadConfigFile_Debug)
{
    ReadConfigFile_Debug();
}
//
//TEST_F(CoreTest,ReadConfigFile_Modules)
//{
//   ASSERT_EQ(c0._modules[0].first->ID,"solar");
//   ASSERT_EQ(c0._modules[1].first->ID,"shadows");
//}