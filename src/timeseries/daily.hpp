#pragma once

#include "timeseries.hpp"
#include <string>
#include <algorithm>

namespace daily
{
    timeseries::iterator start_of_day(timeseries& ts, timeseries::iterator& now);
    timeseries::iterator end_of_day(timeseries& ts, timeseries::iterator& now);
    double max(timeseries& ts,timeseries::iterator& now, std::string variable );
    double min(timeseries& ts,timeseries::iterator& now, std::string variable);
}