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
#ifndef FILTERRATE_H
#define FILTERRATE_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterRate
 * @ingroup processing
 * @author Thomas Egger - Mathias Bavay
 * @date   2011-04-19
 * @brief Rate of change filter.
 * Calculate the change rate (ie: slope) between two points, if it is above a user given value, reject the point.
 *  - If one argument is provided, it is interpreted as the absolute value of the maximum permissible rate of change (per seconds). This means that
 *    every point where <em>|local_rate_of_change| \> argument</em> is rejected
 *  - If two arguments are provided, they are interpreted as the minimum and the maximum (respectively) permissible rate of change (per seconds). This means that
 *    every point where <em>local_rate_of_change \< argument1 AND local_rate_of_change \> argument2</em> is rejected
 *
 * @code
 * TA::filter1	= rate
 * TA::arg1	= -0.01 0.015
 * @endcode
 */
class FilterRate : public FilterBlock {
	public:
		FilterRate(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector<std::string>& vec_args);
		double min_rate_of_change, max_rate_of_change;
};

} //end namespace

#endif
