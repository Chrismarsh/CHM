#pragma once
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_linalg.h>
#include <boost/throw_exception.hpp>
#include <exception.hpp>
#include <iostream>
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
     * Allocating A,b,x, and p every time this is called adds up.
     * So the sz can be pre-set to pre allocate all the arrays.
     * However, this assumes, with no checking, that in calls to operator()
     * sample_points.size() == sz
     */
    thin_plate_spline(size_t sz);

    /**
    * Spline of the sample_points at the query_point location.
    * \param sample_points Tuple of x,y,z values that comprise the sample points from which to interpolate from
    * \param query_point Tuple of x,y,z value that is the point to interpolate to
    * \return Interpolated value at the query_point
    */
    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point);
private:
    double *A ;
    double pi;
    double c; //euler constant
    double *B; // known values - constant value of 0 goes in b[size-1]
    double weight;
    bool static_size;
    size_t size;
    gsl_matrix_view m;
    gsl_vector_view b;
    gsl_vector* x;
    gsl_permutation* p;
};
