#pragma once

#include "interp_alg_base.hpp"
class spline : public interp_alg_base
{
public:
    spline();
    ~spline();
    double operator()(station_list const&  stations, mesh_elem& elem, std::string variable, boost::shared_ptr<interp_visitor> visitor);
       
};
