#pragma once

#include "interpolation.hpp"
class spline : public interpolation_base
{
public:
    spline();
    ~spline();
    double operator()(station_list const&  stations, mesh_elem& elem, std::string variable);
       
};
