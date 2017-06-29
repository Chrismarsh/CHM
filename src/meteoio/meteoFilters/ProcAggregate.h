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
#ifndef PROCAGGREGATE_H
#define PROCAGGREGATE_H

#include <meteoio/meteoFilters/WindowedFilter.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcAggregate
 * @ingroup processing
 * @brief Data aggregation.
 * This aggregates the input data over the defined window with the defined aggregation algorithm. The aggregation
 * algorithm must be declared first and can be any of the following:
 *    + min: return the minimum value of the whole window;
 *    + max: return the maximum value of the whole window;
 *    + mean: return the mean of the whole window;
 *    + median: return the median of the whole window;
 *    + wind_avg: Wind vector averaging. CURRENTLY, THIS FILTER DOES NOT WORK PROPERLY (the first parameter is correctly calculated but the second one uses the modified output of the first one and therefore is WRONG).
 * 
 * Remarks:
 * - nodata values are excluded from the aggregation
 * - Two other arguments are expected (both have to be fullfilled for the filter to start operating):
 *   - minimal number of points in window
 *   - minimal time interval spanning the window (in seconds)
 * - the two arguments may be preceded by the keywords "left", "center" or "right", indicating the window position
 * - the keyword "soft" maybe added (this is highly recommended), if the window position is allowed to be adjusted to the data present
 * @code
 *          VW::filter3 = AGGREGATE
 *          VW::arg3    = MEAN soft left 4 14400 ;(14400 seconds time span for the left leaning window)
 * @endcode
 */

class ProcAggregate : public WindowedFilter {
	public:
		ProcAggregate(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		typedef enum AGGREGATE_TYPE {
			min_agg,
			max_agg,
			mean_agg,
			median_agg,
			wind_avg_agg
		} aggregate_type;
		
		void parse_args(std::vector<std::string> vec_args);
		double calc_min(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const;
		double calc_max(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const;
		double calc_mean(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const;
		double calc_median(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const;
		double calc_wind_avg(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const;
		
		aggregate_type type;
};

} //end namespace

#endif
