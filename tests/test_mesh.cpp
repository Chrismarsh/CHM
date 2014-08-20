

#include "../src/mesh/triangulation.h"
#include "gtest/gtest.h"

class MeshTest : public testing::Test
{
protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
            
      
    }


};

TEST_F(MeshTest,ThrowsOnInvalidFile)
{
    triangulation t0;
    EXPECT_ANY_THROW(t0.from_file("invalid_file.txt"));
}
