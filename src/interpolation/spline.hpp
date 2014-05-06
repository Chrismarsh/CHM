#pragma once

#include "interp_alg_base.hpp"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
class spline : public interp_alg_base
{
public:
    spline();
    ~spline();
    double operator()(station_list&  stations, mesh_elem& elem,boost::shared_ptr<interp_visitor> visitor, boost::shared_ptr<global> global_param);
       
};
