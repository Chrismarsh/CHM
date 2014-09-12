#pragma once

#include "interp_base.hpp"
#include <armadillo>

/**
* \class inv_dist
* Inverse distance weighting interpolation. Requires a minimum of two sets of sample_points
*/
class inv_dist : public interp_base
{
public:
    inv_dist();
    ~inv_dist();

    /**
    * IDW of the sample_points at the query_point location.
    * \param sample_points Tuple of x,y,z values that comprise the sample points from which to interpolate from
    * \param query_point Tuple of x,y,z value that is the point to interpolate to
    * \return Interpolated value at the query_point
    */
    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point);
           
};