#pragma once


#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>

#include <tbb/concurrent_vector.h>

#include <math.h>
#include <iostream>
#include <string>

#include "interpolation_base.hpp"

#include "station.hpp"
#include "triangle.h"

#include "inv_dist.hpp"
#include "spline.hpp"


typedef tbb::concurrent_vector< boost::shared_ptr<station> > station_list;

class interpolation
{
public:
    interpolation();
    ~interpolation();
    double operator()(std::string method, station_list const& stations, mesh_elem& elem, std::string variable);
};





