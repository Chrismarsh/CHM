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
#ifndef __FILTERMAD_H__
#define __FILTERMAD_H__

#include <meteoio/meteoFilters/WindowedFilter.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterMAD
 * @ingroup processing
 * @date   2011-11-01
 * @brief Median Absolute Deviation.
 * Values outside of median ± 3 σ_MAD are rejected. The σ_MAD is calculated as follows:\n
 * <center>\f$ \sigma_{MAD} = K \cdot \mathop{median_i} \left( \left| X_i - \mathop{median_j} ( X_j ) \right| \right) \f$ with \f$ K = \Phi^{-1}; \Phi = 0.6745 \f$ </center>\n
 * See http://en.wikipedia.org/wiki/Median_absolute_deviation for more information about the Mean Absolute Deviation.
 * @code
 * Valid examples for the io.ini file:
 *          TA::filter1 = mad
 *          TA::arg1    = soft left 1 1800  (1800 seconds time span for the left leaning window)
 *          RH::filter1 = mad
 *          RH::arg1    = 10 600            (strictly centered window spanning 600 seconds and at least 10 points)
 * @endcode
 */

class FilterMAD : public WindowedFilter {
	public:
		FilterMAD(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void MAD_filter_point(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end, double &value);
		void parse_args(std::vector<std::string> vec_args);
};

} //end namespace

#endif
