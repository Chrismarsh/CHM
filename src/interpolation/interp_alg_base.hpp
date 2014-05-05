#pragma once

#include <tbb/concurrent_vector.h>
#include <string>

#include "station.hpp"
#include "triangle.h"
#include "logger.h"
#include "exception.hpp"
#include "global.hpp"



typedef tbb::concurrent_vector< boost::shared_ptr<station> > station_list;

class interp_visitor
{
public:
       interp_visitor(){};
       virtual ~interp_visitor(){};
       
       //is run as a modifier to the data pre interpolation
       virtual double lower(mesh_elem & m,  boost::shared_ptr<station>  s, boost::shared_ptr<global> global_param)=0;
       
       //is run as a modifier to the data post interpolation
       virtual double raise(double value, mesh_elem & m, boost::shared_ptr<global> global_param)=0;

};

class interp_alg_base
{
public:
       interp_alg_base(){};
       virtual ~interp_alg_base(){};
       virtual double operator()(station_list&  stations, mesh_elem& elem, boost::shared_ptr<interp_visitor> visitor, boost::shared_ptr<global> global_param){return -9999;};
};




