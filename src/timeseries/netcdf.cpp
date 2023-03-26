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
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Variable not initialized: " + var));
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
void netcdf::open_GEM(const std::string &file)
{
    _data.open(file.c_str(), netCDF::NcFile::read);

    // gem netcdf files have 1 coordinate, datetime

//    for (auto& itr : _data.getVars())
//    {
//        LOG_DEBUG << itr.first;
//    }
//
//    LOG_DEBUG << "-----";
//    for (auto& itr : _data.getCoordVars())
//    {
//        LOG_DEBUG << itr.first;
//    }
//    LOG_DEBUG << "-----";
//    for (auto& itr : _data.getDims())
//    {
//        LOG_DEBUG << itr.first;
//    }
    auto coord_vars = _data.getCoordVars();

    // a few NC have time as a variable and not a coordinate variable so look for time/datetime there
    if(coord_vars.size() == 0)
    {
        CHM_THROW_EXCEPTION(forcing_error,"Netcdf file does not have a coordinate variable defined.");
    }
    _datetime_field = "datetime";

    try
    {
        _datetime_length = coord_vars[_datetime_field].getDim(_datetime_field).getSize();
    }
    catch (netCDF::exceptions::NcNullGrp& e)
    {
        _datetime_field = "time";
        try {
            _datetime_length = coord_vars[_datetime_field].getDim(_datetime_field).getSize();
        }
        catch (netCDF::exceptions::NcNullGrp& e)
        {
            LOG_ERROR << "Tried datetime and time, coord not found";
            throw e;
        }

    }


    // if we don't have at least two timesteps, we can't figure out the model internal timestep length (dt)
    if(_datetime_length == 1)
    {
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("There needs to be at least 2 timesteps in order to determine model dt."));
    }
//
//    if(coord_vars.size() > 1)
//    {
//        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Too many coordinate variables."));
//    }
//
//    for(auto itr: _data.getCoordVars())
//    {
//        _datetime_field = itr.first;
//        _datetime_length = itr.second.getDim(_datetime_field).getSize();
//    }

    netCDF::NcVar times = _data.getVar(_datetime_field);

    if(times.getType().getName() != "int64")
    {
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Datetime dimension not in int64 format"));
    }


    //load in the time offsets
    int64_t* dt = new int64_t[_datetime_length];
    times.getVar(dt);

    //figure out what the epoch is
    // we are expecting the units attribute data to look like
    // hours since 2018-01-05 01:00:00
    std::string epoch;
    auto a = times.getAtt("units");
    a.getValues(epoch);


    if( epoch.find("hours") != std::string::npos )
    {
        LOG_DEBUG << "Found epoch offset = hours";
        _timestep = boost::posix_time::hours(1);
        dt_unit = "hours";
    }
    else if( epoch.find("minutes") != std::string::npos )
    {
        LOG_DEBUG << "Found epoch offset = minutes";
        _timestep = boost::posix_time::minutes(1);
        dt_unit = "minutes";
    }
    else if(( epoch.find("seconds") != std::string::npos ))
    {
        LOG_DEBUG << "Found epoch offset = seconds";
        _timestep = boost::posix_time::seconds(1);
        dt_unit = "seconds";

    } else
    {
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Unknown datetime epoch offset unit."));
    }

    std::vector<std::string> strs;
    boost::split(strs, epoch, boost::is_any_of(" "));

    if(strs.size() != 4)
    {

        //might be in iso format (2017-08-13T01:00:00)
        if(strs.size() != 3)
        {
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Epoch did not split properly, unknown units/ Epoch as read was: " + epoch));
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

    //need to handle a start that is different from our epoch
    // e.g., the epoch might be 'hours since 2021-01-01 00:00:00',
    // but timestep 1 is "5 hours" making the start 2021-01-01 05:00:00
    _start = _epoch + _timestep * dt[0];

    //figure out what the end of the timeseries is
    _end = _epoch + _timestep * dt[_datetime_length-1];


    //get our dt, assuming constant dt throughout the nc file
    _timestep *= dt[1]-dt[0];

    // go through all the timesteps and ensure a consistent timesteping
    // best to spend the time up front for this check than to get 90% into a sim and have it die
    size_t pred_timestep = dt[0];
    for(size_t i=1;  // intentional
         i<_datetime_length; i++)
    {
        pred_timestep += (dt[1]-dt[0]);

        if( dt[i] != pred_timestep)
        {
            std::stringstream expected;
            expected << _epoch + _timestep * pred_timestep;
            std::stringstream got;
            got << _epoch + _timestep * dt[i];


            CHM_THROW_EXCEPTION(forcing_error, "The timesteps in the netcdf file are not constant. At timestep " +
                                                   std::to_string(i) + " offset " + std::to_string(pred_timestep) + " was expected but found " +
                                std::to_string(dt[i]) + ".\n Expected=" + expected.str() + "\n Got=" + got .str()
                                );

        }
    }

    LOG_DEBUG << "NetCDF epoch is " << _epoch;
    LOG_DEBUG << "NetCDF start is " << _start;
    LOG_DEBUG << "NetCDF end is " << _end;
    LOG_DEBUG << "NetCDF timestep is " << _timestep;

    auto dims = _data.getDims();
    for(auto itr : dims)
    {
        if(itr.first == "xgrid_0")
            xgrid = itr.second.getSize();
        else if(itr.first == "ygrid_0")
            ygrid = itr.second.getSize();
    }

    LOG_DEBUG << "NetCDF grid is " << xgrid << " (x) by " << ygrid << " (y)";


}

size_t netcdf::get_ntimesteps()
{
    return _datetime_length;
}

boost::posix_time::time_duration netcdf::get_dt()
{
    return _timestep;
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
    return get_var("HGT_P0_L1_GST",0);
}

std::set<std::string> netcdf::get_variable_names()
{
    if(_variable_names.empty())
    {
        auto vars = _data.getVars();

        std::vector<std::string> exclude = {"datetime","leadtime", "reftime", "HGT_P0_L1_GST", "gridlat_0", "gridlon_0", "xgrid_0", "ygrid_0"};

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

    netcdf::data array(boost::extents[ygrid][xgrid]);
    auto itr = vars.find(var);
    itr->second.getVar(startp,countp, array.data());

    double fill_value = get_fillvalue(itr->second);

    for(size_t i =0; i< array.shape()[0]; i++)
    {
        for(size_t j =0; j< array.shape()[1]; j++)
        {
            if (array[i][j] == fill_value)
                array[i][j] = std::nan("nan");
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

netcdf::data netcdf::get_lat()
{
    return get_var2D("gridlat_0");
}
netcdf::data netcdf::get_lon()
{
    return get_var2D("gridlon_0");
}

size_t netcdf::get_xsize()
{
    return xgrid;
}
size_t netcdf::get_ysize()
{
    return ygrid;
}

double netcdf::get_lat(size_t x, size_t y)
{
    return get_var2D("gridlat_0",x,y);
}
double netcdf::get_lon(size_t x, size_t y)
{
    return get_var2D("gridlon_0",x,y);
}

double netcdf::get_z(size_t x, size_t y)
{
    return get_var("HGT_P0_L1_GST", 0, x, y);
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


    netcdf::data array(boost::extents[ygrid][xgrid]);

    auto itr = vars.find(var);
    itr->second.getVar(startp,countp, array.data());

    double fill_value = get_fillvalue(itr->second);


    for(size_t i =0; i< array.shape()[0]; i++)
    {
        for(size_t j =0; j< array.shape()[1]; j++)
        {
            if (array[i][j] == fill_value)
                array[i][j] = std::nan("nan");
        }
    }

    return array;

}

double netcdf::get_var(std::string var, boost::posix_time::ptime timestep, size_t x, size_t y)
{
    auto diff = timestep - _start; // a duration

    auto offset = diff.total_seconds() / _timestep.total_seconds();

    return get_var(var, offset,x,y);
}
netcdf::data netcdf::get_var(std::string var, boost::posix_time::ptime timestep)
{
    auto diff = timestep - _start; // a duration

    auto offset = diff.total_seconds() / _timestep.total_seconds();

    return get_var(var, offset);
}