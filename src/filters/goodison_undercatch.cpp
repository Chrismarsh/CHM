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

#include "goodison_undercatch.hpp"
REGISTER_FILTER_CPP(goodison_undercatch);

goodison_undercatch::goodison_undercatch(config_file cfg)
  : filter_base("goodison_undercatch", cfg)
{

}
goodison_undercatch::~goodison_undercatch()
{

}
void goodison_undercatch::init()
{
    //look at the config data to determine what we are modifying
    precip_var = cfg.get<std::string>("precip_var");
    wind_var = cfg.get<std::string>("wind_var");
}
void goodison_undercatch::process(std::shared_ptr<station>& station)
{

    double data = (*station)[precip_var];
    double u = (*station)[wind_var];
    //trap missing data, just ignore it.
    if(data != 0) //  CE * 0 will still be zero, but if u is NaN, we will NaN our precip, which we don't want if p = 0
    {
        if( !is_nan(data) && !is_nan(u))
        {
            double CR = 100.00 - 0.44*u*u-1.98*u; // in %
            CR /= 100.0; // fraction
            data /= CR; // pg 46 S4.9.3
        } else
        {
            data = -9999;
        }
    }

    (*station)[precip_var]=data;

}
