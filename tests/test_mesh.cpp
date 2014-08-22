

#include "../src/mesh/triangulation.h"
#include "gtest/gtest.h"

class MeshTest : public testing::Test
{
protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
        EXPECT_NO_THROW(t.from_file("tin_30mdem_30mtol_nodes.csv"));
      
    }

    triangulation t;
};

TEST_F(MeshTest,ThrowsOnInvalidFile)
{
    triangulation t0;
    EXPECT_ANY_THROW(t0.from_file("invalid_file.txt"));
}

TEST_F(MeshTest,ToFileVTU)
{
    EXPECT_NO_THROW(t.to_vtu("test.vtu"));
}