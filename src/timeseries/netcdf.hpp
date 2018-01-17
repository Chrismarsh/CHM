#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp> // for boost::posix
#include <netcdf>
#include <string>

#include "logger.hpp"
#include "exception.hpp"


// handles loading in a netcdf file and loading on the fly nc data into the timeseries class
class netcdf
{

public:
    typedef std::vector< boost::posix_time::ptime > date_vec;

    netcdf();
    ~netcdf();

    boost::posix_time::time_duration get_dt(); //get timestep
    boost::posix_time::ptime get_start();
    boost::posix_time::ptime get_end();


    void open(const std::string& file);
private:

    netCDF::NcFile _data; // main netcdf file

    std::string _datetime_field; // name of the datetime field, the unlimited dimension
    size_t _datetime_length; //number of records

    boost::posix_time::ptime _start, _end;
    date_vec _datetime; //holds all the datetimes

    boost::posix_time::time_duration _timestep; // we will multiply this later to get proper offset

};