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
#ifndef __WINDOWEDFILTER_H__
#define __WINDOWEDFILTER_H__

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  WindowedFilter
 * @brief
 * @author Thomas Egger
 * @date   2011-01-22
 */

class WindowedFilter : public FilterBlock {
	public:
		enum Centering {
			left,   ///< left centered window
			center, ///< centered window
			right   ///< right centered window
		};

		WindowedFilter(const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec) = 0;

	protected:
		static unsigned int get_centering(std::vector<std::string>& vec_args);
		const std::vector<const MeteoData*>& get_window(const size_t& index,
		                                                const std::vector<MeteoData>& ivec);
		bool get_window_specs(const size_t& index, const std::vector<MeteoData>& ivec,
		                      size_t& start, size_t& end) const;

		Duration min_time_span;
		Centering centering;
		size_t last_start, last_end;
		size_t min_data_points;

	private:
		std::vector<const MeteoData*> vec_window;

	protected:
		bool is_soft; //placed here for alignement
};

} //end namespace

#endif
