#pragma once
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <boost/throw_exception.hpp>
#include <exception.hpp>

#include "interp_base.hpp"

/**
* \class thin_plate_spline
*
* Thin plate spline with tensions interpolation
*/
class thin_plate_spline : public interp_base
{
public:
    thin_plate_spline();
    ~thin_plate_spline();

    /**
    * Spline of the sample_points at the query_point location.
    * \param sample_points Tuple of x,y,z values that comprise the sample points from which to interpolate from
    * \param query_point Tuple of x,y,z value that is the point to interpolate to
    * \return Interpolated value at the query_point
    */
    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point);

};
