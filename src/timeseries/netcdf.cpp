

#include "netcdf.hpp"

netcdf::netcdf()
{

}
netcdf::~netcdf()
{

}
void netcdf::open(const std::string& file)
{
    _data.open(file.c_str(), netCDF::NcFile::read);

    //gem netcdf files have 1 coordinate, datetime
    auto coord_vars = _data.getCoordVars();

    if(coord_vars.size() > 1)
    {

    }
    for(auto itr: _data.getCoordVars())
    {
        _datetime_field = itr.first;
        _datetime_length = itr.second.getDim(_datetime_field).getSize();
    }

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
    }
    else if( epoch.find("minutes") != std::string::npos )
    {
        LOG_DEBUG << "Found epoch offset = minutes";
        _timestep = boost::posix_time::minutes(1);
    }
    else if(( epoch.find("seconds") != std::string::npos ))
    {
        LOG_DEBUG << "Found epoch offset = seconds";
        _timestep = boost::posix_time::seconds(1);

    } else
    {
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Unknown datetime epoch offset unit."));
    }

    std::vector<std::string> strs;
    boost::split(strs, epoch, boost::is_any_of(" "));

    if(strs.size() != 4)
    {
        BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Epoch did not split into 4 sections, unknown units/ Epoch as read was: " + epoch));
    }

    _start = boost::posix_time::time_from_string(strs[2]+" "+strs[3]);
    LOG_DEBUG << "NetCDF epoch is " << _start;

    //get our dt, assuming constant dt throughout the nc file
    _timestep *= dt[1]-dt[0];

    _end = _start + _timestep*_datetime_length;

    LOG_DEBUG << "NetCDF end is " << _end;
    LOG_DEBUG << "NetCDF timestep is " << _timestep;

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