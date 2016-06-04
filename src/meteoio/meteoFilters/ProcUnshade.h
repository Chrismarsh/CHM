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
#ifndef __PROCUNSHADE_H__
#define __PROCUNSHADE_H__

#include <meteoio/meteoFilters/WindowedFilter.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <vector>
#include <string>
#include <algorithm>

namespace mio {

/**
 * @class  ProcUnshade
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2011-02-07
 * @brief Correct short wave measurements for shading of the sensor.
 *
 * @code
 *
 * @endcode
 */

class ProcUnshade : public WindowedFilter {
	public:
		ProcUnshade(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		bool filterAlbedo(const std::vector<double>& julian, std::vector<double> &data) const;
		bool linInterpolate(const std::vector<double>& julian, std::vector<double> &data) const;
		void parse_args(std::vector<std::string> vec_args);
		static bool MADFilter(const size_t& start_idx, std::vector<double> &vecWindow, std::vector<double> &data);
		static void interpolFill(const size_t& start_idx, const size_t& end_idx, const std::vector<double>& julian, std::vector<double> &data);

		double max_gap; ///<largest albedo gap that will still be interpolated, in days
};

} //end namespace

#endif
