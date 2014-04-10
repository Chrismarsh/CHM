#pragma once

#include "interpolation.hpp"


class inv_dist : public interpolation_base
{
public:
    inv_dist();
    ~inv_dist();
    double operator()(station_list const&  stations, mesh_elem& elem, std::string variable);
           
};