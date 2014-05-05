#pragma once

#include "interp_alg_base.hpp"
#include <armadillo>

class inv_dist : public interp_alg_base
{
public:
    inv_dist();
    ~inv_dist();
    double operator()(station_list&  stations, mesh_elem& elem, boost::shared_ptr<interp_visitor> visitor, boost::shared_ptr<global> global_param);
           
};