#pragma once
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <cmath>

#include "exception.hpp"

class interp_base
{
public:
    //sample_points is a vector of xyz points
    //query_xy are the x,y coords of the point to interp to
    virtual double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_points )
    {
        return 0.0;
    };

    virtual ~interp_base(){};
    interp_base(){};

};