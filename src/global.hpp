#pragma once
//#include "core.h"
#include <boost/date_time/posix_time/posix_time.hpp>
//#include <boost/move/unique_ptr.hpp>
#include <string>

#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

#include "interpolation.h"
#include "station.hpp"
/**
 * Basin wide parameters such as transmissivity, solar elevation, solar aspect, etc.
 * 
 **/
class global
{
    //want to let core modify date time, etc without showing a public interface.
    //This is because global gets passed to all modules and a rogue module could do something dumb
    //const doesn't save us as we actually do want to modify things
    friend class core;
    
private:
    //solar elevation in radians
    double _solar_el;
    //solar azimuth in radians
    double _solar_az;
    
    //approximate lat and lon of the basin for solar radiation calculations
    double _lat;
    double _lon;
    
    //approximate UTC offset
    int _utc_offset;
    
    boost::posix_time::ptime _current_date;
    
    void solar_el_az();
    
    //updates the internal variables 
    void update();
    
    var _variables;

protected:
    int _dt; //seconds

public:

    global();
    int year();
    int day();
    interp_alg interp_algorithm;
    /*
     * Montboosth on [1,12]
     */
    int month();
    int hour();
    int min();
    int sec();
    int dt();
    boost::posix_time::ptime posix_time();
    uint64_t posix_time_int();
    double solar_el();
    double solar_az();
    
    bool first_time_step;
    std::string get_variable(std::string variable);

    //each station where observations are
    std::vector< boost::shared_ptr<station> > stations; //tbb::concurrent_vector

//    template<class T>
//    T get_parameter_value(std::string key);

    pt::ptree parameters;
    
};
//
//template<class T>
//T global::get_parameter_value(std::string key)
//{
//    return parameters.get<T>(key);
//}