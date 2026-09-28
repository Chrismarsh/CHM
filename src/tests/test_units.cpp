#include "gtest/gtest.h"
#include "PhysConst.h"

using namespace PhysConst::units;

/*
 *  Construct units of each type
 *  
 *  Add them together, substract them, no impl mulitply
 *
 *  > < == <= >= += -=
 */ 
template<class A, class B>
concept EqualityComparable =
    requires (A a, B b) {
        { a == b } -> std::convertible_to<bool>;
    };

template<class A, class B>
concept ThreeWayComparable =
    requires (A a, B b) {
        { a <=> b };
    };


class UnitsTest : public ::testing::Test
{
protected:
    struct Example : public base<double,Example> {};
    class Other : public base<double,Other> {};

    
    static_assert(!EqualityComparable<Example,Other>);
    static_assert(!ThreeWayComparable<Example,Other>);
};

TEST_F(UnitsTest,ConstructEmptyIsZero)
{
    Example example;

    EXPECT_EQ(example.value,0.0);
};

TEST_F(UnitsTest,ConstructInitializerBracesFromDouble)
{
    constexpr double init = 4.0;
    Example exam{init};

    EXPECT_EQ(exam.value,init);
};

TEST_F(UnitsTest,PlusEqual)
{
    Example e1{5.0}, e2{6.5};
    
    e1 += e2;

    EXPECT_EQ(e1.value,11.5);
};

TEST_F(UnitsTest,MinusEqual)
{
    Example e1{9.0}, e2{2.5};

    e1 -= e2;

    EXPECT_EQ(e1.value,6.5);
};

TEST_F(UnitsTest,GreaterThan)
{
    Example e1{10.0}, e2{5.5};

    EXPECT_TRUE(e1 > e2);
    EXPECT_TRUE(e1 >= e2);
    EXPECT_TRUE(e1 >= e1);
};

TEST_F(UnitsTest,LessThan)
{
    Example e1{1000.0}, e2{450.0};

    EXPECT_TRUE(e1 < e2);
    EXPECT_TRUE(e1 <= e2);
    EXPECT_TRUE(e1 <= e1);

};

TEST_F(UnitsTest,Equal)
{
    Example e1{3486.35};

    EXPECT_TRUE(e1 == e1);
};
