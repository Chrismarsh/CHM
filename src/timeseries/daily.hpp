#pragma once

#include "timeseries.hpp"
#include <string>
#include <algorithm>

namespace daily
{
    double mean(time_series::iterator& start, time_series::iterator& end, std::string variable );
    double max(time_series::iterator& start, time_series::iterator& end, std::string variable );
    
}