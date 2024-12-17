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
#include <algorithm>
#include <array>

#include "logger.hpp"
#include "exception.hpp"

// handles loading in a netcdf file and loading on the fly nc data into the timeseries class
class netcdf
{

public:
  typedef boost::multi_array<double,2> _ma_data;
  typedef boost::multi_array<double,1> _ma_vec;
  typedef std::shared_ptr<_ma_data> data;
  typedef std::shared_ptr<_ma_vec> vec;
  typedef std::vector< boost::posix_time::ptime > date_vec;

  netcdf();
  ~netcdf();

  /**
   * Returns the model integration timestep
   * @return
   */
boost::posix_time::time_duration get_dt();
  /**
   * Global start time
   * @return
   */
  boost::posix_time::ptime get_start();
  boost::posix_time::ptime get_end();

  std::set<std::string> get_variable_names();
  std::set<std::string> get_coordinate_names();
  void open_GEM(const std::string &file);
  void open(const std::string &file);

  void create(const std::string& file);
  size_t get_xsize() const;
  size_t get_ysize() const;
  size_t get_ntimesteps();
  //returns the lat grid
  vec get_lat();
  double get_lat(size_t x, size_t y);

  /**
   * 2D lat grid, such as the rotated HRPS grids
   * @return
   */
  data get_lat2D();
  //returns the 1D lon grid
  vec get_lon();
 /**
 * 2D lat grid, such as the rotated HRPS grids
 * @return
 */
  data get_lon2D();
  double get_lon(size_t x, size_t y);

  //gets z information
  data get_z();
  double get_z(size_t x, size_t y);

  /**
  * Get the CF _FillValue if present, defaults to -9999.0 if that attribute is not present.
  * Currently assumes double!
  * @param var
  * @return
  */
  double get_fillvalue(const netCDF::NcVar& var);

  data get_var(std::string var, size_t timestep);
  data get_var(std::string var, boost::posix_time::ptime timestep);
  double get_var(std::string var, size_t timestep, size_t x, size_t y);
  double get_var(std::string var, boost::posix_time::ptime timestep, size_t x, size_t y);

  /**
   * Finds a coordinate by a standard_name, e.g., "time"
   * @param standard_name The standar_name attr to search for
   * @return
   */
  std::string find_coord_by_standard_name(const std::string& standard_name) const;

  /**
 * Finds a dimension by a standard_name, e.g., "latitude"
 * @param standard_name The standar_name attr to search for
 * @return
 */
  std::string find_dim_by_standard_name(const std::string& standard_name) const;

  /**
   * Finds a variable by looking in various attributes. Defaults to {"standard_name", "long_name"}
   * @param search String to search for. Exact match
   * @param attrs_to_search Vector of attributes to search
   * @return
   */
 std::string find_var_by_attr(const std::string& search, const std::vector<std::string>& attrs_to_search={"standard_name", "long_name"}) const;

  void add_dim1D(const std::string& var, size_t length);
  void create_variable1D(const std::string& var,  size_t length);
  void put_var1D(const std::string& var, size_t index, double value);
  /**
   * Some data, such as lat/long do not have a time component are only 2D data. This allows loading those data.
   * @param var
   * @return
   */
  data get_var2D(std::string var);
  double get_var1D(std::string var, size_t index);
  double get_var2D(std::string var, size_t x, size_t y);
  std::string get_var_standard_name(const std::string& variable) const;

 // do we have a z data variable?
 bool missing_z();

  /**
   * Returns the dimensionality of the spatial coordinates. 1D or 2D
   * @return
   */
int get_coord_dimensionality() const;
  netCDF::NcFile& get_ncfile();
private:

  netCDF::NcFile _data; // main netcdf file
  std::string _datetime_field; // name of the datetime field, the unlimited dimension
  std::string _lat_field, _lon_field; //name of lat and long fields
  std::string _elevation_field; // name of elevation field
  size_t xgrid, ygrid;

 //dimension of the spatial coordinate. 1 or 2 is currently supported
  int _spatial_coord_dim;

  std::set<std::string> _variable_names; //set of variables this nc file provides

  size_t _datetime_length; //number of records

  boost::posix_time::ptime _start, _epoch, _end;

  boost::posix_time::time_duration _epoch_offset_unit; // we will multiply this later to get proper offset
  boost::posix_time::time_duration _delta_t;  // model timestep size. whole integer values, in seconds. e.g, 3600s
  bool _is_open;

  // if we are missing z coords that need to be estimated later
  bool _missing_z;


  //if we are creating variables
  std::vector<netCDF::NcDim> _dimVector; //we need this dimension var to create new variables

};