#include "math/Bisection.hpp"
#include "gtest/gtest.h"
#include <cmath>
#include <functional>

using namespace math;

class BisectionTest : public ::testing::Test
{

protected:
    static constexpr auto foo = [](const double x) {
        return x * x - 4.0;
    };
};

TEST_F(BisectionTest,NoRootGuaranteed)
{
    auto result = bisection(foo,-1.0,1.0);

    EXPECT_TRUE(std::isnan(result.value));
    EXPECT_EQ(result.root,Root::NoRootGuaranteed);
};

TEST_F(BisectionTest,DoesNotConverge)
{
    constexpr size_t matched_exponent = 3;
    static constexpr auto slow = [](const double x) {
        return std::log(x) - matched_exponent;
    };

    auto result = bisection(slow,1.0,200.0,1e-12,10);

    EXPECT_EQ(result.root,Root::DidNotConverge);
    EXPECT_EQ(result.iterations,10);

    result = bisection(slow,std::exp(matched_exponent)*0.8,std::exp(matched_exponent)*1.2,1e-3);

    EXPECT_EQ(result.root,Root::Found);

    EXPECT_EQ(result.value,std::exp(matched_exponent)); 
}

TEST_F(BisectionTest,XSquaredRoot)
{

    auto result = bisection(foo,-1.0,3.0);

    EXPECT_EQ(result.root,Root::Found);
    EXPECT_NEAR(result.value,2.0,1e-10);
};
