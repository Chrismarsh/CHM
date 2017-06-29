/***********************************************************************************/
/*  Copyright 2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of MeteoIO.
    MeteoIO is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MeteoIO is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with MeteoIO.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef FILTERNOCHANGE_H
#define FILTERNOCHANGE_H

#include <meteoio/meteoFilters/WindowedFilter.h>

#include <vector>
#include <string>

namespace mio {

/**
 * @class FilterNoChange
 * @ingroup processing
 * @brief This filter removes periods showing insufficient changes. 
 * It searches for time periods in which the value of the certain variable doesn't change by looking at the variance. 
 * It expects the following arguments, in order to define the data window where the variance will be computed: 
 * - the optional keyword "soft"
 * - the optional window centering option (left, center or right)
 * - the minimal number of points in window
 * - the minimal time interval spanning the window (in seconds)
 * 
 * For example:
 * @code
 *          HS::filter1 = NO_CHANGE
 *          HS::arg1    = soft left 1 1800 (1800 seconds time span for the left leaning window)
 * 
 *          TA::filter1 = NO_CHANGE
 *          TA::arg1    = 10 600          (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 * 
 * @author Anna-Maria Tilg
 * @date   2015-12-04
 */

class FilterNoChange : public WindowedFilter {
	public:
		FilterNoChange(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
};

} //end namespace

#endif
