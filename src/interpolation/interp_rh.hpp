#pragma once
#include <string>
#include <boost/shared_ptr.hpp>
#include "interp_alg_base.hpp"
#include "interp_2d.hpp"
#include "variable_map.hpp"


class interp_rh 
{
public:
    interp_rh();
    ~interp_rh();
    void operator()(std::string method, mesh_elem& m, station_list& stations, boost::shared_ptr<global> global_param);
};

class LLRA_rh_var: public interp_visitor
{
public:
    double lower(mesh_elem& m,  boost::shared_ptr<station>  s, boost::shared_ptr<global> global_param);
    double raise(double value, mesh_elem& m, boost::shared_ptr<global> global_param);
private:
    double get_lambda_rate(int month);
};


