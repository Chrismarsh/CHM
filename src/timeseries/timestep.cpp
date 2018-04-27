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


#include "timestep.hpp"


timestep::timestep(const boost::shared_ptr<timestep> src)
{
   
    _itrs = itr_map(src->_itrs);

    _date_itr = date_variable::iterator(src->_date_itr);
}

timestep::timestep()
{
#if USE_SPARSEHASH
   _itrs.set_empty_key("");
#endif
}

timestep::~timestep()
{

}

std::string timestep::to_string()
{
    std::stringstream s;
    s << get_posix() << "\t";
    for (timestep::itr_map::const_iterator itr = _itrs.begin(); itr != _itrs.end(); itr++)
    {
        s << boost::lexical_cast<std::string>(*(itr->second)) << std::string("\t");
    }

    return s.str();
}

int timestep::month()
{
    return _date_itr->date().month();
}


int timestep::day()
{
    return _date_itr->date().day();

}

int timestep::year()
{
    return _date_itr->date().year();
}

int timestep::hour()
{
    return boost::posix_time::to_tm(*_date_itr).tm_hour;
}
int timestep::min()
{
    return boost::posix_time::to_tm(*_date_itr).tm_min;
}
int timestep::sec()
{
    return boost::posix_time::to_tm(*_date_itr).tm_sec;
}


boost::gregorian::date timestep::get_gregorian()
{
    boost::gregorian::date date;

    date = boost::gregorian::from_string(boost::lexical_cast<std::string>(_date_itr->date()));
    return date;

}

boost::posix_time::ptime timestep::get_posix()
{
    return *_date_itr;
}

bool timestep::has(const std::string &variable)
{
    auto res = _itrs.find(variable);

    if(res == _itrs.end())
        return false;
    else
        return true;
}

double timestep::get(const std::string &variable)
{
#ifdef SAFE_CHECKS
    auto res = _itrs.find(variable);
    if(res == _itrs.end())
        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + variable + " does not exist."));
#endif

    return _itrs[variable][0];

}

void timestep::set(const std::string &variable, const double &value)
{
    //this is slower than just getting the value, but it prevents subtle bugs where
    // a module tries to create a variable it didn't allocate in a provides call
#ifdef SAFE_CHECKS
    auto res = _itrs.find(variable);
    if(res == _itrs.end())
        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + variable + " does not exist."));
#endif

    _itrs[variable][0] = value;

}

timestep::variable_vec::iterator timestep::get_itr(const std::string &varName)
{
    auto res = _itrs.find(varName);
    if(res == _itrs.end())
        BOOST_THROW_EXCEPTION( forcing_lookup_error() << errstr_info("Variable " + varName + " does not exist."));

    return res->second;

}

