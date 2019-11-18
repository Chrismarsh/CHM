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

#include <stddef.h>
#include "station.hpp"

station::~station()
{

}

station::station()
{
    _x = 0;
    _y = 0;
    _z = 0.0;

}

station::station(std::string ID, double x, double y, double elevation, std::set<std::string> variables )
{
    _ID = ID;
    _x = x;
    _y = y;
    _z = elevation;
    init(variables);
}

double& station::operator[](const uint64_t& hash)
{
    return _timestep_data[hash];
}
double& station::operator[](const std::string& variable)
{
    return _timestep_data[variable];
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
   return strm << std::fixed << "ID=" << s._ID << "; (x,y,z)=(" << s._x << "," << s._y << "," << s._z << "); ";
}

void station::init(std::set<std::string> variables)
{
    _timestep_data.init(variables);
}

boost::gregorian::date station::get_gregorian()
{
    boost::gregorian::date date;

    date = boost::gregorian::from_string(boost::lexical_cast<std::string>(_current_ts.date()));
    return date;

}
boost::posix_time::ptime station::get_posix()
{
    return _current_ts;
}

void station::set_posix(boost::posix_time::ptime ts )
{
    _current_ts = ts;
}

bool station::has(const std::string &variable)
{
    return false;
//    auto res = _itrs.find(variable);
//
//    if(res == _itrs.end())
//        return false;
//    else
//        return true;
}


int station::month()
{
    return _current_ts.date().month();
}


int station::day()
{
    return _current_ts.date().day();

}

int station::year()
{
    return _current_ts.date().year();
}

int station::hour()
{
    return boost::posix_time::to_tm(_current_ts).tm_hour;
}
int station::min()
{
    return boost::posix_time::to_tm(_current_ts).tm_min;
}
int station::sec()
{
    return boost::posix_time::to_tm(_current_ts).tm_sec;
}
