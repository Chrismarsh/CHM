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


#include "interpolation.hpp"
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

TEST_F(InterpTest,spline2)
{
    thin_plate_spline s;
    std::vector<boost::tuple<double,double,double> > xy;

    xy.push_back( boost::make_tuple(-1276639.4142831599,1408220.6433826166,22.241299818717572));
    xy.push_back( boost::make_tuple(-1276628.96002623, 1408213.5776356135, 22.423794697169313));
    xy.push_back( boost::make_tuple(-1276628.8896492834,1408225.6645281466,22.301020204404736));

    auto query = boost::make_tuple(-1276633.6294519969,1408220.6575855566,2306.0533040364585);

    double z0 =  s(xy,query);
    std::cout << z0 << std::endl;

    double result = s(xy,query);

    ASSERT_DOUBLE_EQ(result,22.301293951848109);


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

TEST_F(InterpTest,interpolation_class_static_size2)
{
    interpolation s(interp_alg::tpspline,3);
    std::vector<boost::tuple<double,double,double> > xy;

    xy.push_back( boost::make_tuple(-1276639.4142831599,1408220.6433826166,22.241299818717572));
    xy.push_back( boost::make_tuple(-1276628.96002623, 1408213.5776356135, 22.423794697169313));
    xy.push_back( boost::make_tuple(-1276628.8896492834,1408225.6645281466,22.301020204404736));

    auto query = boost::make_tuple(-1276633.6294519969,1408220.6575855566,2306.0533040364585);

    double z0 =  s(xy,query);
    std::cout << z0 << std::endl;

    double result = s(xy,query);

    ASSERT_DOUBLE_EQ(result,22.301293951848109);


}
