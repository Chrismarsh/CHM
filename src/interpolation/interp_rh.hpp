#pragma once
#include <string>
#include <boost/shared_ptr.hpp>
#include "interp_alg_base.hpp"
#include "interp_2d.hpp"
class interp_rh 
{
public:
    interp_rh();
    ~interp_rh();
    void operator()(std::string method, mesh_elem& m, station_list& stations, std::string variable);
};

class LLRA_rh_var: public interp_visitor
{
public:
    double lower(mesh_elem& m, std::string temperature_id, boost::shared_ptr<station>  s);
    double raise(double value, mesh_elem& m, std::string temperature_id);
private:
    double get_lambda_rate(int month);
    int _month;
};


