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



#include "netcdf.hpp"

netcdf::netcdf()
{
    _is_open = false;
    _datetime_field="";
}
netcdf::~netcdf()
{

}
 void netcdf::add_dim1D(const std::string& var, size_t length)
 {
     auto nTri = _data.addDim(var, length);
     _dimVector.push_back(nTri);
 }
void netcdf::create_variable1D( const std::string& var, size_t length)
{
    //only create the dim and variables once
    try
    {
        add_dim1D("tri_id",length);

    }
    catch(netCDF::exceptions::NcNameInUse& e)
    {

    }

    try
    {
        auto nc_var = _data.addVar(var.c_str(), netCDF::ncDouble, _dimVector);
    }
    catch(netCDF::exceptions::NcNameInUse& e)
    {

    }


}

netCDF::NcFile& netcdf::get_ncfile()
{
    return _data;
}

void netcdf::put_var1D(const std::string& var, size_t index, double value)
{
    auto vars = _data.getVars();

    auto itr = vars.find(var);

    std::vector<size_t> startp,countp;
    startp.push_back(index);
    countp.push_back(1);

    try
    {
        itr->second.putVar(startp,countp,&value);
    }
    catch(netCDF::exceptions::NcBadId& e)
    {
        CHM_THROW_EXCEPTION(forcing_error, "Variable not initialized: " + var);
    }



}

void netcdf::create(const std::string& file)
{
    _data.open(file.c_str(), netCDF::NcFile::replace);

}
void netcdf::open(const std::string &file)
{
    _data.open(file.c_str(), netCDF::NcFile::read);
}

std::string netcdf::find_coord_by_standard_name(const std::string& standard_name) const
{
    for (auto& itr : _data.getCoordVars())
    {
        auto var = _data.getVar(itr.first);
        auto standardNameAtt = var.getAtt("standard_name");
        if(standardNameAtt.isNull())
            continue;

        std::string std_name;
        standardNameAtt.getValues(std_name);
        if(std_name == standard_name)
        {
            return itr.first;
        }
    }

    // if we get here, we didn't find what we were looking for
    CHM_THROW_EXCEPTION(forcing_error, "Could not find coordinate variable for " + standard_name);

}

std::string netcdf::find_dim_by_standard_name(const std::string& standard_name) const
{
    for (auto& itr : _data.getDims())
    {
        auto var = _data.getVar(itr.first);
        netCDF::NcVarAtt standardNameAtt;
        try
        {
            standardNameAtt = var.getAtt("standard_name");

            if(standardNameAtt.isNull())
                continue;
        }
        catch(netCDF::exceptions::NcException& e)
        {
            // this doesn't even have a standard_name attr
            continue;
        }

        std::string std_name;
        standardNameAtt.getValues(std_name);
        if(std_name == standard_name)
        {
            return itr.first;
        }
    }

    // if we get here, we didn't find what we were looking for
    CHM_THROW_EXCEPTION(forcing_error, "Could not find coordinate variable for " + standard_name);

}

void netcdf::open_GEM(const std::string &file)
{
    _data.open(file.c_str(), netCDF::NcFile::read);

    // gem netcdf files have 1 coordinate, datetime

    SPDLOG_DEBUG("NC getVars:");
    for (auto& itr : _data.getVars())
    {
        SPDLOG_DEBUG("\t"+itr.first);
    }

    SPDLOG_DEBUG("NC getCoordVars");
    for (auto& itr : _data.getCoordVars())
    {
        SPDLOG_DEBUG("\t"+itr.first);
    }
    SPDLOG_DEBUG("Nc getDims");
    for (auto& itr : _data.getDims())
    {
        SPDLOG_DEBUG("\t"+itr.first);
    }

    auto coord_vars = _data.getCoordVars();

    // a few NC have time as a variable and not a coordinate variable so look for time/datetime there
    if(coord_vars.size() == 0)
    {
        CHM_THROW_EXCEPTION(forcing_error,"Netcdf file does not have a coordinate variable defined.");
    }

    _datetime_field = find_coord_by_standard_name("time");
    SPDLOG_DEBUG("Found time coordinate variable: " + _datetime_field);

    try
    {
        _datetime_length = coord_vars[_datetime_field].getDim(_datetime_field).getSize();
    }
    catch (netCDF::exceptions::NcNullGrp& e)
    {
        SPDLOG_ERROR("Could not find datetime coordinate");
        throw;
    }

//
//    if(coord_vars.size() > 1)
//    {
//        CHM_THROW_EXCEPTION(forcing_error, "Too many coordinate variables.");
//    }
//
//    for(auto itr: _data.getCoordVars())
//    {
//        _datetime_field = itr.first;
//        _datetime_length = itr.second.getDim(_datetime_field).getSize();
//    }

    netCDF::NcVar times = _data.getVar(_datetime_field);

    //load in the time offsets
    auto* dt = new int64_t[_datetime_length];

    if(times.getType().getName() == "int64")
    {
        times.getVar(dt);
    }
    else  if(times.getType().getName() == "double")
    {
        auto* tmp_dt = new double[_datetime_length];
        times.getVar(tmp_dt);

        for(int i = 0; i < _datetime_length; i++)
        {
            dt[i] = static_cast<int64_t>(tmp_dt[i]);
        }

        delete[] tmp_dt;
    }
    else
    {
        CHM_THROW_EXCEPTION(forcing_error, "Datetime dimension not in int64 or double type. Type is: " + times.getType().getName());
    }

    //figure out what the epoch is
    // we are expecting the units attribute data to look like
    // hours since 2018-01-05 01:00:00
    std::string epoch;
    auto a = times.getAtt("units");

    try
    {
        // SPDLOG_DEBUG(times.getAtt("units").getType().getName());
        a.getValues(epoch);
    }
    catch(...)
    {
        SPDLOG_ERROR("Datetime attributes are likely string attributes. Did you write this netcdf with xarray with engine=h5netcdf? Try engine=netcdf4");
        CHM_THROW_EXCEPTION(forcing_error, "Datetime string attributes not supported");
    }



    if( epoch.find("hours") != std::string::npos )
    {
        SPDLOG_DEBUG("Found epoch offset = hours");
        _epoch_offset_unit = boost::posix_time::hours(1);
    }
    if( epoch.find("days") != std::string::npos )
    {
        SPDLOG_DEBUG("Found epoch offset = days");
        _epoch_offset_unit = boost::posix_time::hours(24);
    }
    else if( epoch.find("minutes") != std::string::npos )
    {
        SPDLOG_DEBUG("Found epoch offset = minutes");
        _epoch_offset_unit = boost::posix_time::minutes(1);
    }
    else if(( epoch.find("seconds") != std::string::npos ))
    {
        SPDLOG_DEBUG("Found epoch offset = seconds");
        _epoch_offset_unit = boost::posix_time::seconds(1);

    } else
    {
        CHM_THROW_EXCEPTION(forcing_error, "Unknown datetime epoch offset unit: " + epoch);
    }

    std::vector<std::string> strs;
    boost::split(strs, epoch, boost::is_any_of(" "));

    if(strs.size() != 4)
    {
        //might be in iso format (2017-08-13T01:00:00)
        if(strs.size() != 3)
        {
            CHM_THROW_EXCEPTION(forcing_error, "Epoch did not split properly, unknown units/ Epoch as read was: " + epoch);
        }

        //If it's 3, means there is a T b/w date and time, remove it.
        std::string s = strs[2];
        auto tpos = s.find("T");
        if (tpos != std::string::npos)
        {
            s.replace(s.find("T"),1," ");
        }

        // midnight times can be reported without the 00:00 suffix. If we get this far and don't have : in the epoch
        // then we need to add it
        tpos = s.find(":");
        if (tpos == std::string::npos)
        {
            s = s + " 00:00:00";
        }

        try
        {
            _epoch = boost::posix_time::time_from_string(s);
        }
        catch(boost::bad_lexical_cast& e)
        {
            CHM_THROW_EXCEPTION(forcing_error, "Unable to parse netcdf epoch time " + s);
        }
    }
    else
    {
        _epoch = boost::posix_time::time_from_string(strs[2]+" "+strs[3]);
    }

    //get our dt, expectation is that it is a constant dt throughout the nc file
    try
    {
        auto attr = coord_vars[_datetime_field];

        size_t dt = 0;
        // has to be from the variable not the coord to get access to the attr
        _data.getVar(_datetime_field).getAtt("delta_t").getValues(&dt);

        std::string dt_units;
        // this getValues is different than the rest
        // https://docs.unidata.ucar.edu/netcdf-cxx/current/classnetCDF_1_1NcAtt.html#a07ba5f59a1d9a1e1d0eca6adf133796c
        _data.getVar(_datetime_field).getAtt("delta_t_units").getValues(dt_units);

        if( dt_units != "s")
        {
            CHM_THROW_EXCEPTION(forcing_error, "time:delta_t_units must be s");
        }

        _delta_t = boost::posix_time::seconds(dt);

    }
    catch(forcing_error& e)
    {
        CHM_THROW_EXCEPTION(forcing_error, e.what());
    }
    catch(...)
    {
        // no attribute, try to compute it manually
        // if we don't have at least two timesteps, we can't figure out the model internal timestep length (dt)
        if(_datetime_length == 1)
        {
            CHM_THROW_EXCEPTION(forcing_error,"There needs to be at least 2 timesteps in order to determine model dt or the time coordinate needs to have the attribute 'delta_t:<step in seconds'.");
        }
        _delta_t = boost::posix_time::seconds(dt[1]-dt[0]);
    }

    //need to handle a start that is different from our epoch
    // e.g., the epoch might be 'hours since 2021-01-01 00:00:00',
    // but timestep 1 is "5 hours" making the start 2021-01-01 05:00:00
    _start = _epoch + _epoch_offset_unit * dt[0];

    //figure out what the end of the timeseries is
    _end = _epoch + _epoch_offset_unit * dt[_datetime_length-1];


    // go through all the timesteps and ensure a consistent timesteping
    // best to spend the time up front for this check than to get 90% into a sim and have it die
    size_t pred_timestep = dt[0];
    for(size_t i=1;  // intentional
         i < _datetime_length; i++)
    {
        pred_timestep =  pred_timestep + _delta_t.total_seconds();

        if( dt[i] != pred_timestep)
        {
            std::stringstream expected;
            expected << _epoch + _epoch_offset_unit * pred_timestep;
            std::stringstream got;
            got << _epoch + _epoch_offset_unit * dt[i];


            CHM_THROW_EXCEPTION(forcing_error, "The timesteps in the netcdf file are not constant. At timestep " +
                                                   std::to_string(i) + " offset " + std::to_string(pred_timestep) + " was expected but found " +
                                std::to_string(dt[i]) + ".\n Expected=" + expected.str() + "\n Got=" + got .str()
                                );

        }
    }

    SPDLOG_DEBUG("NetCDF epoch is {}", boost::posix_time::to_simple_string(_epoch));
    SPDLOG_DEBUG("NetCDF start is {}", boost::posix_time::to_simple_string(_start));
    SPDLOG_DEBUG("NetCDF end is {}",boost::posix_time::to_simple_string(_end));
    SPDLOG_DEBUG("NetCDF timestep is {}", boost::posix_time::to_simple_string(_delta_t));

    // CF convention assumes that dim and coord have the same name,
    _lat_field = find_dim_by_standard_name("latitude");
    _lon_field = find_dim_by_standard_name("longitude");
    _elevation_field = find_var_by_attr("geopotential_height");

    xgrid = _data.getDim(_lon_field).getSize();
    ygrid = _data.getDim(_lat_field).getSize();

    SPDLOG_DEBUG("NetCDF grid is {} (x) by {} (y)", xgrid, ygrid);

    auto lat_dim = _data.getVar(find_coord_by_standard_name("latitude")).getDimCount();
    auto lon_dim = _data.getVar(find_coord_by_standard_name("longitude")).getDimCount();

    if(lat_dim != lon_dim)
    {
        CHM_THROW_EXCEPTION(forcing_error, "Latitude and longitude dimensionality do not match");
    }

    _spatial_coord_dim = lat_dim;
    SPDLOG_DEBUG("Coord dimensionality is {}D",_spatial_coord_dim);

    if(_spatial_coord_dim > 2)
    {
        CHM_THROW_EXCEPTION(forcing_error, "Latitude and longitude dimensionality exceeds 2D");
    }


}

size_t netcdf::get_ntimesteps()
{
    return _datetime_length;
}

boost::posix_time::time_duration netcdf::get_dt()
{
    return _delta_t;
}

boost::posix_time::ptime netcdf::get_start()
{
    return _start;
}
boost::posix_time::ptime netcdf::get_end()
{
    return _end;
}

netcdf::data netcdf::get_z()
{
    return get_var(_elevation_field, 0);
}


std::string netcdf::find_var_by_attr(const std::string& search, const std::vector<std::string>& attrs_to_search) const
{
    // even though we only return 1 value, we need to make sure we didn't match multiple
    std::vector<std::string> result;

    if(attrs_to_search.empty())
    {
        CHM_THROW_EXCEPTION(forcing_error, "Empty attribute list provided");
    }

    for (auto& [fst, snd] : _data.getVars())
    {
        auto var = _data.getVar(fst);

        // check all requested attributes
        for(const auto& jtr:attrs_to_search)
        {
            netCDF::NcVarAtt search_att;
            try
            {
                // no standard_name, ignore it
                search_att = var.getAtt(jtr);
                if(search_att.isNull())
                    continue;
            }
            catch(netCDF::exceptions::NcException& e)
            {
                // no standard_name, ignore it
                continue;
            }


            std::string std_name;
            search_att.getValues(std_name);
            if(std_name == search)
            {
                result.push_back(fst);
            }
        }
    }

    if (result.empty())
    {
        // we didn't find what we were looking for
        CHM_THROW_EXCEPTION(forcing_error, "Could not find variable by searching attrs for " + search);
    }
    else if(result.size() > 1)
    {
        CHM_THROW_EXCEPTION(forcing_error, "Multiple variables found for attr " + search);
    }

    return result.at(0);


}

std::set<std::string> netcdf::get_variable_names()
{
    if(_variable_names.empty())
    {
        auto vars = _data.getVars();

        std::vector<std::string> exclude = {"datetime","leadtime", "reftime", "HGT_P0_L1_GST", _lat_field, _lon_field};

        for (auto itr: vars)
        {
            auto v = itr.first;
            //don't return the above variables as they are geo-spatial vars
            if (std::find(exclude.begin(), exclude.end(), v) == exclude.end())
            {
                _variable_names.insert(v);
            }
        }
    }

    return _variable_names;
}

std::set<std::string> netcdf::get_coordinate_names()
{
    std::set<std::string> names;
    auto vars = _data.getCoordVars();

    for (auto itr: vars)
    {
        auto v = itr.first;
        names.insert(v);
    }

    return names;

}

double netcdf::get_fillvalue(const netCDF::NcVar& var)
{
    double fill_value=-9999;
    try {
        auto fill_value_attr  = var.getAtt("_FillValue");
        fill_value_attr.getValues(&fill_value);
    }catch(netCDF::exceptions::NcException& e)
    {
        // no CF _FillValue, use -9999 default
    }
    return fill_value;
}

double netcdf::get_var1D(std::string var, size_t index)
{
    std::vector<size_t> startp, countp;

    startp.push_back(index);
    countp.push_back(1);


    auto vars = _data.getVars();

    auto itr = vars.find(var);
    double data=-9999.0;
    itr->second.getVar(startp,countp,&data);


    double fill_value = get_fillvalue(itr->second);

    if( data == fill_value)
    {
        data = std::nan("nan");
    }

    return data;
}

netcdf::data netcdf::get_var2D(std::string var)
{
    std::vector<size_t> startp, countp;

    startp.push_back(0);
    startp.push_back(0);

    countp.push_back(ygrid);
    countp.push_back(xgrid);

    auto vars = _data.getVars();

    netcdf::data array = std::make_shared<_ma_data>(_ma_data(boost::extents[ygrid][xgrid]));
    auto itr = vars.find(var);
    itr->second.getVar(startp,countp, array->data());

    double fill_value = get_fillvalue(itr->second);

    for(size_t i =0; i< array->shape()[0]; i++)
    {
        for(size_t j =0; j< array->shape()[1]; j++)
        {
            if ((*array)[i][j] == fill_value)
                (*array)[i][j] = std::nan("nan");
        }
    }

    return array;
}

double netcdf::get_var2D(std::string var, size_t x, size_t y)
{
    std::vector<size_t> startp, countp;

    startp.push_back(y);
    startp.push_back(x);

    countp.push_back(1);
    countp.push_back(1);

    auto vars = _data.getVars();

    double val=-9999;

    auto itr = vars.find(var);
    itr->second.getVar(startp,countp, &val);

    double fill_value = get_fillvalue(itr->second);

    if(val == fill_value)
        val = std::nan("nan");

    return val;
}

std::string netcdf::get_var_standard_name(const std::string& variable) const
{
    auto var = _data.getVar(variable);
    auto standardNameAtt = var.getAtt("standard_name");
    if(standardNameAtt.isNull())
        return "";

    std::string std_name;
    standardNameAtt.getValues(std_name);

    return std_name;
}

netcdf::data netcdf::get_lat2D()
{
    return get_var2D(_lat_field);
}

netcdf::data netcdf::get_lon2D()
{
    return get_var2D(_lon_field);
}

netcdf::vec netcdf::get_lat()
{
    auto vars = _data.getVars();
    auto itr = vars.find(_lat_field);
    netcdf::vec array = std::make_shared<_ma_vec>(_ma_vec(boost::extents[ygrid]));
    itr->second.getVar({0},{ygrid}, array->data());

    double fill_value = get_fillvalue(itr->second);

    std::transform(array->begin(), array->end(), array->begin(),
                   [fill_value](double val)
                   {
                       return (val == fill_value) ? std::nan("") : val;
                   });


    return array;
}
netcdf::vec netcdf::get_lon()
{
    auto vars = _data.getVars();
    auto itr = vars.find(_lon_field);
    netcdf::vec array = std::make_shared<_ma_vec>(_ma_vec(boost::extents[xgrid]));
    itr->second.getVar({0},{xgrid}, array->data());

    double fill_value = get_fillvalue(itr->second);

    std::transform(array->begin(), array->end(), array->begin(),
                   [fill_value](double val)
                   {
                       return (val == fill_value) ? std::nan("") : val;
                   });


    return array;
}

size_t netcdf::get_xsize() const
{
    return xgrid;
}
size_t netcdf::get_ysize() const
{
    return ygrid;
}

double netcdf::get_lat(size_t x, size_t y)
{
    return get_var2D(_lat_field,x,y);
}
double netcdf::get_lon(size_t x, size_t y)
{
    return get_var2D(_lon_field,x,y);
}

double netcdf::get_z(size_t x, size_t y)
{
    return get_var(_elevation_field, 0, x, y);
}

int netcdf::get_coord_dimensionality() const
{
    return _spatial_coord_dim;
}


double netcdf::get_var(std::string var, size_t timestep, size_t x, size_t y)
{
    std::vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(y);
    startp.push_back(x);

    countp.push_back(1);
    countp.push_back(1);
    countp.push_back(1);

    // Read the data one record at a time.
    startp[0] = timestep;

    auto vars = _data.getVars();

    double val=-9999;

    auto itr = vars.find(var);
#pragma omp critical
    {
        itr->second.getVar(startp, countp, &val);
    }
    double fill_value = get_fillvalue(itr->second);

    if(val == fill_value)
        val = std::nan("nan");

    return val;
}

netcdf::data netcdf::get_var(std::string var, size_t timestep)
{
    std::vector<size_t> startp, countp;
    startp.push_back(0);
    startp.push_back(0);
    startp.push_back(0);

    countp.push_back(1);
    countp.push_back(ygrid);
    countp.push_back(xgrid);

    // Read the data one record at a time.
    startp[0] = timestep;

    auto vars = _data.getVars();


    netcdf::data array = std::make_shared<_ma_data>(boost::extents[ygrid][xgrid]);

    auto itr = vars.find(var);
    itr->second.getVar(startp,countp, array->data());

    double fill_value = get_fillvalue(itr->second);

    for(size_t i =0; i< array->shape()[0]; i++)
    {
        for(size_t j =0; j< array->shape()[1]; j++)
        {
            if ((*array)[i][j] == fill_value)
                (*array)[i][j] = std::nan("nan");
        }
    }

    return array;

}

double netcdf::get_var(std::string var, boost::posix_time::ptime timestep, size_t x, size_t y)
{
    auto diff = timestep - _start; // a duration

    auto offset = diff.total_seconds() / _epoch_offset_unit.total_seconds();

    return get_var(var, offset,x,y);
}
netcdf::data netcdf::get_var(std::string var, boost::posix_time::ptime timestep)
{
    auto diff = timestep - _start; // a duration

    auto offset = diff.total_seconds() / _epoch_offset_unit.total_seconds();

    return get_var(var, offset);
}