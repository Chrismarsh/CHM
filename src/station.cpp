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
    _face = 0;

}

station::station(std::string ID, double x, double y, double elevation)
{
    _ID = ID;
    _x = x;
    _y = y;
    _z = elevation;
    _obs = new timeseries();
    _itr = _obs->begin();
    _face = 0;

}

size_t station::closest_face()
{
    return _face;
}
void station::set_closest_face(size_t face)
{
    _face = face;
}
void station::open(std::string file)
{
    try
    {
        _obs = new timeseries();
        _obs->open(file);

        _itr = _obs->begin();
    }
    catch (exception_base &e)
    {
        //e << errstr_info( std::string("Station:") + _ID ); //hope like hell we've got the ID at this point
        throw;
    }
}


timeseries* station::raw_timeseries()
{
    return _obs;
}

timeseries::date_vec station::date_timeseries()
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

double station::x()
{
    return _x;
}

double station::y()
{
    return _y;
}

void station::x(double x)
{
    _x = x;
}

void station::y(double y)
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
void station::tofile(std::string file)
{
    _obs->to_file(file);
}