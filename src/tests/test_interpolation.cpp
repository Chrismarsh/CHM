
#include "interpolation.h"
#include "logger.hpp"
#include <vector>
#include <boost/tuple/tuple.hpp>

#include "gtest/gtest.h"
#include <stdlib.h>
#include <string>
#include <utility>


class InterpTest : public testing::Test
{

protected:

    virtual void SetUp()
    {
        logging::core::get()->set_logging_enabled(false);
    }

};

TEST_F(InterpTest,spline)
{
    thin_plate_spline s;
    std::vector<boost::tuple<double,double,double> > xy;

    xy.push_back( boost::make_tuple(69.,76.,20.820));
    xy.push_back( boost::make_tuple(59.,64.,10.910 ));
    xy.push_back( boost::make_tuple(75.,52.,10.380 ));
    xy.push_back( boost::make_tuple(86.,73.,14.600 ));
    xy.push_back( boost::make_tuple(88.,53.,10.560 ));

    auto query = boost::make_tuple(69.,67.,0.);

    double z0 =  s(xy,query);
    std::cout << z0 << std::endl;

    double result = fabs(s(xy,query) - 15.795);

    ASSERT_LT(result,1);

}

TEST_F(InterpTest,spline_forcedsize)
{
    thin_plate_spline s(5);
    std::vector<boost::tuple<double,double,double> > xy;

    xy.push_back( boost::make_tuple(69.,76.,20.820));
    xy.push_back( boost::make_tuple(59.,64.,10.910 ));
    xy.push_back( boost::make_tuple(75.,52.,10.380 ));
    xy.push_back( boost::make_tuple(86.,73.,14.600 ));
    xy.push_back( boost::make_tuple(88.,53.,10.560 ));

    auto query = boost::make_tuple(69.,67.,0.);

    std::cout << s(xy,query) << std::endl;
    double result = fabs(s(xy,query) - 15.795);

    ASSERT_LT(result,1);

}

TEST_F(InterpTest,interpolation_class)
{
    interpolation s(interp_alg::tpspline);
    std::vector<boost::tuple<double,double,double> > xy;

    xy.push_back( boost::make_tuple(69.,76.,20.820));
    xy.push_back( boost::make_tuple(59.,64.,10.910 ));
    xy.push_back( boost::make_tuple(75.,52.,10.380 ));
    xy.push_back( boost::make_tuple(86.,73.,14.600 ));
    xy.push_back( boost::make_tuple(88.,53.,10.560 ));

    auto query = boost::make_tuple(69.,67.,0.);

    std::cout << s(xy,query) << std::endl;
    double result = fabs(s(xy,query) - 15.795);

    ASSERT_LT(result,1);

}
TEST_F(InterpTest,interpolation_class_static_size)
{
    interpolation s(interp_alg::tpspline,5);
    std::vector<boost::tuple<double,double,double> > xy;

    xy.push_back( boost::make_tuple(69.,76.,20.820));
    xy.push_back( boost::make_tuple(59.,64.,10.910 ));
    xy.push_back( boost::make_tuple(75.,52.,10.380 ));
    xy.push_back( boost::make_tuple(86.,73.,14.600 ));
    xy.push_back( boost::make_tuple(88.,53.,10.560 ));

    auto query = boost::make_tuple(69.,67.,0.);

    std::cout << s(xy,query) << std::endl;
    double result = fabs(s(xy,query) - 15.795);

    ASSERT_LT(result,1);

}