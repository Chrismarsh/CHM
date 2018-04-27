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

#include "variable_map.hpp"
#include "logger.hpp"
#include "exception.hpp"
var::var()
{
    
}
void var::init_from_file(std::string path)
{
    LOG_DEBUG << "Reading variable from file: " << path;
    
    //tmp to load it up until file read code is implimented
    var_hashmap::accessor a;
    _varmap.insert(a,"RH"); a->second = "rh";
    _varmap.insert(a,"Tair"); a->second = "t";
    _varmap.insert(a,"timestep"); a->second = "datetime";
}

var::~var()
{
    
}

std::string var::operator()(std::string variable)
{
    var_hashmap::accessor a;
    if(!_varmap.find(a, variable))
    {
        BOOST_THROW_EXCEPTION(forcing_error()
                                << errstr_info("Unable to find " + variable));
    }   
    return a->second;
}
