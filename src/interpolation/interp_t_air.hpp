#pragma once
#include <string>
#include <boost/shared_ptr.hpp>
#include "interp_alg_base.hpp"
#include "interp_2d.hpp"
class interp_t_air 
{
public:
    interp_t_air();
    ~interp_t_air();
    void operator()(std::string method, mesh_elem& m, station_list& stations, std::string variable);
};

class LLRA_const : public interp_visitor
{
public:
    double lower(mesh_elem& m, std::string temperature_id, boost::shared_ptr<station>  s);
    double raise(double value, mesh_elem& m, std::string temperature_id);
};

class LLRA_var : public interp_visitor
{
public:
    double lower(mesh_elem& m, std::string temperature_id, boost::shared_ptr<station>  s);
    double raise(double value, mesh_elem& m, std::string temperature_id);
private:
    double get_lapse_rate(int month);
    int _month;
};

