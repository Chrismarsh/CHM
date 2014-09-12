#pragma once

#include "timeseries.hpp"
#include <string>
#include <algorithm>

namespace daily
{
    double mean(timeseries::iterator& start, timeseries::iterator& end, std::string variable );
    double max(timeseries::iterator& start, timeseries::iterator& end, std::string variable );
    
}