//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

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
    typedef boost::multi_array<double,1> vec;
    typedef std::vector< boost::posix_time::ptime > date_vec;

    netcdf();
    ~netcdf();

    boost::posix_time::time_duration get_dt(); //get timestep
    boost::posix_time::ptime get_start();
    boost::posix_time::ptime get_end();

    std::set<std::string> get_variable_names();
    void open_GEM(const std::string &file);
    void open(const std::string &file);

    void create(const std::string& file);
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


    data get_var(std::string var, size_t timestep);
    data get_var(std::string var, boost::posix_time::ptime timestep);
    double get_var(std::string var, size_t timestep, size_t x, size_t y);
    double get_var(std::string var, boost::posix_time::ptime timestep, size_t x, size_t y);

    void add_dim1D(const std::string& var, size_t length);
    void create_variable1D(const std::string& var,  size_t length);
    void put_var1D(const std::string& var, size_t index, double value);
    /**
     * Some data, such as lat/long do not have a time component are only 2 data. This allows loading those data.
     * @param var
     * @return
     */
    data get_var2D(std::string var);
    double get_var1D(std::string var, size_t index);

    double get_var2D(std::string var, size_t x, size_t y);

    date_vec get_datevec();


    netCDF::NcFile& get_ncfile();
private:

    netCDF::NcFile _data; // main netcdf file
    std::string _datetime_field; // name of the datetime field, the unlimited dimension
    std::string _lat, _lon; //name of lat and long fields
    size_t xgrid, ygrid;

    std::set<std::string> _variable_names; //set of variables this nc file provides

    size_t _datetime_length; //number of records

    boost::posix_time::ptime _start, _epoch, _end;
    date_vec _datetime; //holds all the datetimes
    std::string dt_unit;

    boost::posix_time::time_duration _timestep; // we will multiply this later to get proper offset

    bool _is_open;


    //if we are creating variables
    std::vector<netCDF::NcDim> _dimVector; //we need this dimension var to create new variables

};