#pragma once
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <boost/throw_exception.hpp>
#include <exception.hpp>

#include "interp_base.hpp"

class thin_plate_spline : public interp_base
{
public:
    thin_plate_spline();
    ~thin_plate_spline();
    double operator()(std::vector< boost::tuple<double,double,double> > sample_points, boost::tuple<double,double,double> query_points );

};
