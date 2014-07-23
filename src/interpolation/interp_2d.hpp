#pragma once


#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>

#include <tbb/concurrent_vector.h>

#include <math.h>
#include <iostream>
#include <string>



#include "station.hpp"
#include "triangulation.h"

#include "interp_alg_base.hpp"
#include "inv_dist.hpp"
#include "spline.hpp"
#include "global.hpp"

typedef tbb::concurrent_vector< boost::shared_ptr<station> > station_list;

class interp_2d
{
public:
    interp_2d();
    ~interp_2d();
    double operator()(std::string method, mesh_elem& elem, station_list&  stations, boost::shared_ptr<interp_visitor> visitor, boost::shared_ptr<global> global_param);
};




