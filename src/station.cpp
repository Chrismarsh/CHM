#include <stddef.h>
#include "station.hpp"

station::~station()
{
    delete _obs;
}

station::station()
{
    _x = 0;
    _y = 0;
    _z = 0.0;
    _obs = NULL;

}

station::station(std::string ID, size_t x, size_t y, double elevation)
{
    _ID = ID;
    _x = x;
    _y = y;
    _z = elevation;
    _obs = new time_series();
    _itr = _obs->begin();

}

void station::open(std::string file)
{
    try
    {
        _obs = new time_series();
        _obs->open(file);

        _itr = _obs->begin();
    }
    catch (exception_base &e)
    {
        //e << errstr_info( std::string("Station:") + _ID ); //hope like hell we've got the ID at this point
        throw;
    }
}

time_series::date_vec station::date_timeseries()
{
    return _obs->get_date_timeseries();
}

size_t station::timeseries_length()
{
    return _obs->get_timeseries_length();
}

std::vector<std::string> station::list_variables()
{
    return _obs->list_variables();
}


timestep& station::now()
{
    return *_itr;
}

double station::get(std::string variable)
{
    return _itr->get(variable);
}

void station::reset_itrs()
{
    _itr = _obs->begin();
}

bool station::next()
{
    ++_itr;
    if (_itr == _obs->end())
        return false;
    else
        return true;
}

size_t station::x()
{
    return _x;
}

size_t station::y()
{
    return _y;
}

void station::x(size_t x)
{
    _x = x;
}

void station::y(size_t y)
{
    _y = y;
}

double station::z()
{
    return _z;
}

void station::z(double elevation)
{
    _z = elevation;
}

void station::ID(std::string ID)
{
    _ID = ID;
}

std::string station::ID()
{
    return _ID;
}

std::ostream &operator<<(std::ostream &strm, const station &s)
{
    if (s._obs)
        return strm << "ID=" << s._ID << " (x,y,z)=(" << s._x << "," << s._y << "," << s._z << ") ," << "forcing=" << s._obs->get_opened_file();
    else
        return strm << "ID=" << s._ID << "; (x,y,z)=(" << s._x << "," << s._y << "," << s._z << "); " << "forcing=Not opened";
}