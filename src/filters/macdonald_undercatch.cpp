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

//
// Created by chris on 18/11/15.
//

#include "macdonald_undercatch.hpp"
REGISTER_FILTER_CPP(macdonald_undercatch);


macdonald_undercatch::macdonald_undercatch(config_file cfg)
  : filter_base("macdonald_undercatch", cfg)
{

}
macdonald_undercatch::~macdonald_undercatch()
{

}
void macdonald_undercatch::init()
{
    //look at the config data to determine what we are modifying
    var = cfg.get<std::string>("variable");
}
void macdonald_undercatch::process(std::shared_ptr<station>& station)
{

    double data = (*station)[var];
    double u = (*station)["u"];
    //trap missing data, just ignore it.
    if( !is_nan(data) && !is_nan(u))
    {
        data /= (1.010 * exp(-0.09*u));
    } else
    {
        data = -9999;
    }

    (*station)[var]=data;

}
