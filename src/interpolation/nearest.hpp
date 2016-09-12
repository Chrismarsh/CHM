#pragma once

#include "interp_base.hpp"
#include <armadillo>

/**
* \class nearest
* Simple interpolation that applies 100% weight to the nearest station. Requires exactly 1 station.
*/
class nearest : public interp_base
{
public:
    nearest();
    ~nearest();

    /**
    * Nearest sample point
    * \param sample_points Tuple of x,y,z values that comprise the sample points from which to interpolate from
    * \param query_point Tuple of x,y,z value that is the point to interpolate to
    * \return Interpolated value at the query_point
    */
    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point);
           
};