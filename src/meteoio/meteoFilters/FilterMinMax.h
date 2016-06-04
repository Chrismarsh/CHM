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
#ifndef __FILTERMINMAX_H__
#define __FILTERMINMAX_H__

//#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterMinMax
 * @ingroup processing
 * @brief Min/Max range filter.
 * @author Thomas Egger
 * @date   2011-01-02
 * Reject all values greater than the max or smaller than the min. Remarks:
 * - two arguments have to be provided, min and max (in SI)
 * - the keyword "soft" maybe added, in such a case all data greater than the max would be assigned
 * the maximum permissible value and all data smaller than the min would be assigned the minimum permissible value
 * or an optional extra set of two user provided values (see example below)
 * @code
 * TA::filter1	= min_max
 * TA::arg1	= 230 330
 * ISWR::filter1	= min_max
 * ISWR::arg1	= soft 8 1500 0 1498
 * @endcode
 *
 */

class FilterMinMax : public FilterBlock {
	public:
		FilterMinMax(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);

		double min_val, max_val;
		double min_soft, max_soft;
		bool is_soft;
};

} //end namespace

#endif
