#pragma once
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <cmath>

#include "exception.hpp"

/**
* \class interp_base
* Base class for all interpolation schemes to inherent from.
*/
class interp_base
{
public:

    /**
    * The interpolation method must implement this method.
    * \param sample_points Tuple of x,y,z values that comprise the sample points from which to interpolate from
    * \param query_point Tuple of x,y,z value that is the point to interpolate to
    * \return Interpolated value at the query_point
    */
    virtual double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)
    {
        return 0.0;
    };

    virtual ~interp_base(){};
    interp_base(){};

};