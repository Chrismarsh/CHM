#pragma once
//#include "core.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <string>

#include "variable_map.hpp"
/**
 * Basin wide paramters such as transmissivity, solar elevation, solar aspec, etc.
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
    
    //aproximae lat and lon of the basin for solar radiation calculations
    double _lat;
    double _lon;
    
    //approximate UTC offset
    int _utc_offset;
    
    boost::posix_time::ptime _current_date;
    
    void solar_el_az();
    
    //updates the internal variables 
    void update();
    
    var _variables;

public:
    
    int year();
    int day();
    int month();
    int hour();
    int min();
    int sec();
    boost::posix_time::ptime posix_time();
    
    double solar_el();
    double solar_az();
    
    
    std::string get_variable(std::string variable);
    
};