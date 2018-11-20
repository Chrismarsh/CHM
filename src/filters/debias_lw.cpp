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

#include "debias_lw.hpp"
REGISTER_FILTER_CPP(debias_lw);

debias_lw::debias_lw(config_file cfg)
  : filter_base("debias_lw", cfg)
{

}

debias_lw::~debias_lw()
{

}

void debias_lw::init(boost::shared_ptr<station>& station)
{
    //look at the config data to determine what we are modifying
    var = cfg.get<std::string>("variable");
    fac = cfg.get<double>("factor");    // LW correction in W/m2

}
void debias_lw::process(boost::shared_ptr<station>& station)
{
    double data = station->now().get(var);
    if(!is_nan(data))
    {
         data = data + fac;
    }
    
    station->now().set(var,data);
}
