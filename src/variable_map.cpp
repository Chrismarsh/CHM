#include "variable_map.hpp"
#include "logger.hpp"
#include "exception.hpp"
var::var()
{
    
}
void var::init_from_file(std::string path)
{
    LOG_DEBUG << "Reading variable from file: " << path;
    
    //tmp to load it up until file read code is implimented
    var_hashmap::accessor a;
    _varmap.insert(a,"RH"); a->second = "rh";
    _varmap.insert(a,"Tair"); a->second = "t";
    _varmap.insert(a,"timestep"); a->second = "datetime";
}

var::~var()
{
    
}

std::string var::operator()(std::string variable)
{
    var_hashmap::accessor a;
    if(!_varmap.find(a, variable))
    {
        BOOST_THROW_EXCEPTION(forcing_error()
                                << errstr_info("Unable to find " + variable));
    }   
    return a->second;
}
