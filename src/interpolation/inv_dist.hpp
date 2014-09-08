#pragma once

#include "interp_base.hpp"
#include <armadillo>

class inv_dist : public interp_base
{
public:
    inv_dist();
    ~inv_dist();
    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_points);
           
};