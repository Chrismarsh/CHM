/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __FILTERWINDAVG_H__
#define __FILTERWINDAVG_H__

#include <meteoio/meteoFilters/WindowedFilter.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterWindAvg
 * @ingroup processing
 * @author Thomas Egger
 * @date   2011-01-24
 * @brief Wind vector averaging.
 * This calculates the vector average over a user given time period. Each wind vector within this period
 * is added and the final sum is normalized by the number of vectors that have been added. Important information:
 * - nodata values are excluded from the mean
 * - Two arguments expected (both have to be fullfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 *   (centered by default)
 * - the keyword "soft" maybe added, if the window position is allowed to be adjusted to the data present
 * Please notice that this filter can currently ONLY be used on VW and DW. It also makes no sens to use it only on one of
 * these two parameters: both should be processed by this filter.
 * CURRENTLY, THIS FILTER DOES NOT WORK PROPERLY (the first parameter is correctly calculated but the second one
 * uses the modified output of the first one and therefore is WRONG).
 * @code
 * Valid examples for the io.ini file:
 *          VW::filter1 = WIND_AVG
 *          VW::arg1    = soft left 1 1800 (1800 seconds time span for the left leaning window)
 *          VW::filter1 = WIND_AVG
 *          VW::arg1    = 10 600          (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */

class FilterWindAvg : public WindowedFilter {
	public:
		FilterWindAvg(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
		double calc_avg(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end);
};

} //end namespace

#endif
