/***********************************************************************************/
/*  Copyright 2011 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __FILTERSTDDEV_H__
#define __FILTERSTDDEV_H__

#include <meteoio/meteoFilters/WindowedFilter.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <vector>
#include <string>
#include <algorithm>

namespace mio {

/**
 * @class  FilterStdDev
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2011-02-07
 * @brief Standard deviation filter.
 * Values outside of mean ± 2 std_dev are rejected.
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = std_dev
 *          TA::arg1    = soft left 1 1800  (1800 seconds time span for the left leaning window)
 *          RH::filter1 = std_dev
 *          RH::arg1    = 10 6000            (strictly centered window spanning 6000 seconds and at least 10 points)
 * @endcode
 */

class FilterStdDev : public WindowedFilter {
	public:
		FilterStdDev(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
		void getStat(const std::vector<MeteoData>& ivec, const unsigned int& param,
		             const size_t& start, const size_t& end, double& stddev, double& mean);
		static const double sigma; ///<How many times the stddev allowed for valid points
};

} //end namespace

#endif
