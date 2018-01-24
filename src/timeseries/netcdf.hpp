#pragma once

#include <boost/multi_array.hpp>
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
    typedef boost::multi_array<double,2> data;
    typedef std::vector< boost::posix_time::ptime > date_vec;

    netcdf();
    ~netcdf();

    boost::posix_time::time_duration get_dt(); //get timestep
    boost::posix_time::ptime get_start();
    boost::posix_time::ptime get_end();

    std::set<std::string> get_variable_names();
    void open(const std::string& file);

    size_t get_xsize();
    size_t get_ysize();
    size_t get_ntimesteps();
    //returns the lat grid
    data get_lat();
    double get_lat(size_t x, size_t y);

    //returns the long grid
    data get_lon();
    double get_lon(size_t x, size_t y);

    //gets z information
    data get_z();
    double get_z(size_t x, size_t y);

    data get_data(std::string var, size_t timestep);
    data get_data(std::string var, boost::posix_time::ptime timestep);
    double get_var(std::string var, size_t timestep, size_t x, size_t y);

    /**
     * Some data, such as lat/long do not have a time component are only 2 data. This allows loading those data.
     * @param var
     * @return
     */
    data get_var2D(std::string var);

    double get_var2D(std::string var, size_t x, size_t y);

    date_vec get_datevec();

private:

    netCDF::NcFile _data; // main netcdf file
    std::string _datetime_field; // name of the datetime field, the unlimited dimension
    std::string _lat, _lon; //name of lat and long fields
    size_t xgrid, ygrid;

    std::set<std::string> _variable_names; //set of variables this nc file provides

    size_t _datetime_length; //number of records

    boost::posix_time::ptime _start, _end;
    date_vec _datetime; //holds all the datetimes
    std::string dt_unit;

    boost::posix_time::time_duration _timestep; // we will multiply this later to get proper offset

    bool _is_open;

};