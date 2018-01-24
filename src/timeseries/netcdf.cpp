

#include "netcdf.hpp"

netcdf::netcdf()
{
    _is_open = false;
}
netcdf::~netcdf()
{

}
void netcdf::open(const std::string& file)
{
    _data.open(file.c_str(), netCDF::NcFile::read);

    //gem netcdf files have 1 coordinate, datetime

    auto coord_vars = _data.getCoordVars();
    _datetime_field = "datetime";
    _datetime_length = coord_vars[_datetime_field].getDim(_datetime_field).getSize();

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
            BOOST_THROW_EXCEPTION(forcing_error() << errstr_info("Epoch did not split properlys, unknown units/ Epoch as read was: " + epoch));
        }

        //If it's 3, means there is a T b/w date and time, remove it.

        std::string s = strs[2];
        s.replace(s.find("T"),1," ");

        _start = boost::posix_time::time_from_string(s);

    }
    else
    {
        _start = boost::posix_time::time_from_string(strs[2]+" "+strs[3]);
    }


    LOG_DEBUG << "NetCDF epoch is " << _start;

    //get our dt, assuming constant dt throughout the nc file
    _timestep *= dt[1]-dt[0];

    _end = _start + _timestep * dt[_datetime_length-1];

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

    //build up the datetime vector so we can easily use it later
    _datetime.resize(_datetime_length);
    for(size_t i=0; i<_datetime_length;i++)
    {
        _datetime[i] = _start + _timestep*i;
    }

}

size_t netcdf::get_ntimesteps()
{
    return _datetime_length;
}

netcdf::date_vec netcdf::get_datevec()
{
    return _datetime;
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
    return get_data("HGT_P0_L1_GST0",0);
}

std::set<std::string> netcdf::get_variable_names()
{
    if(_variable_names.empty())
    {
        auto vars = _data.getVars();

        std::vector<std::string> exclude = {"HGT_P0_L1_GST0", "gridlat_0", "gridlon_0", "xgrid_0", "ygrid_0"};

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

netcdf::data netcdf::get_var2D(std::string var)
{
    std::vector<size_t> startp, countp;

    startp.push_back(1);
    startp.push_back(1);

    countp.push_back(ygrid);
    countp.push_back(xgrid);

    auto vars = _data.getVars();

    netcdf::data array(boost::extents[ygrid][xgrid]);
    auto itr = vars.find(var);
    itr->second.getVar(startp,countp, array.data());

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
    return get_var("HGT_P0_L1_GST0", 0, x, y);
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
    itr->second.getVar(startp,countp, &val);

    return val;
}

netcdf::data netcdf::get_data(std::string var, size_t timestep)
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


     return array;
}

netcdf::data netcdf::get_data(std::string var, boost::posix_time::ptime timestep)
{
    auto diff = timestep - _start; // a duration

    auto offset = diff.total_seconds() / _timestep.total_seconds();

    return get_data(var, offset);
}