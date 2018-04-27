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


#include "daily.hpp"

namespace daily
{

    double min(timeseries& ts,timeseries::iterator& now, std::string variable)
    {
        auto start = start_of_day(ts,now);
        auto end = end_of_day(ts, now);

        return ts.range_min(start, end, variable);
    }

    double max(timeseries& ts,timeseries::iterator& now, std::string variable)
    {
        auto start = start_of_day(ts,now);
        auto end = end_of_day(ts, now);

        return ts.range_max(start, end, variable);
    }

    timeseries::iterator start_of_day(timeseries& ts, timeseries::iterator& now)
    {

        //today but midnight
        boost::posix_time::ptime time ( boost::gregorian::date(now->year(),now->month(),now->day()),
                                        boost::posix_time::time_duration(0,0,0));

        return ts.find(time);

    }

    timeseries::iterator end_of_day(timeseries& ts, timeseries::iterator& now)
    {

        //today but almost midnight
        boost::posix_time::ptime time ( boost::gregorian::date(now->year(),now->month(),now->day()),
                boost::posix_time::time_duration(23,0,0));

        timeseries::iterator itr = ts.find(time);

        int day = now->day();
        while(itr != ts.end() && itr->day() != day ) //stop at the last element in this day
            ++itr;

        return itr;

    }
}
