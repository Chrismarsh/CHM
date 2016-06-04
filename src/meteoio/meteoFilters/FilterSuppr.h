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
#ifndef __FILTERSUPPR_H__
#define __FILTERSUPPR_H__

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterSuppr
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2013-12-06
 * @brief Suppression filter.
 * Normally, this filter simply reject all values. This is convenient to quickly turn a parameter off
 * without modifying the original data. 
 * But it is also possible to suppress a given fraction of the data at random by providing
 * such fraction as an argument. For example, <i>0.5</i> would ensure that at least <i>50%</i> of the
 * data set contains <i>nodata</i> for this parameter.
 * @code
 * ILWR::filter1     = suppr
 * 
 * TA::filter1       = suppr
 * TA::arg1          = 0.5
 * @endcode
 */

class FilterSuppr : public FilterBlock {
	public:
		FilterSuppr(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);
		
		double range;
};

} //end namespace

#endif
