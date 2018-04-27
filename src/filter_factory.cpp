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

#include "filter_factory.h"


filter_base* filter_factory::get(std::string ID, pt::ptree config)
{
    LOG_VERBOSE << "Filter ID=" << ID;

    filter_base* filter = nullptr;

    if (ID == "macdonald_undercatch")
        filter = new macdonald_undercatch();
    else if (ID == "scale_wind_speed")
        filter = new scale_wind_speed();

    if(filter == nullptr)
    {
        BOOST_THROW_EXCEPTION(module_not_found()
                              << errstr_info( std::string("Filter not found ") + ID)
        );
    }

    filter->ID = ID;
    filter->cfg = config;

    return filter;
}