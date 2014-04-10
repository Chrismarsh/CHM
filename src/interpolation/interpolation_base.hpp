#pragma once

#include <tbb/concurrent_vector.h>
#include <string>

#include "station.hpp"
#include "triangle.h"

typedef tbb::concurrent_vector< boost::shared_ptr<station> > station_list;

class interpolation_base
{
public:
       interpolation_base(){};
       virtual ~interpolation_base(){};
       virtual double operator()(station_list const&  stations, mesh_elem& elem, std::string variable)=0;
};
