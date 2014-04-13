#pragma once

#include "interp_alg_base.hpp"


class inv_dist : public interp_alg_base
{
public:
    inv_dist();
    ~inv_dist();
    double operator()(station_list&  stations, mesh_elem& elem, std::string variable, boost::shared_ptr<interp_visitor> visitor);
           
};