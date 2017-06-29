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
#ifndef FILTERMAX_H
#define FILTERMAX_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterMax
 * @ingroup processing
 * @author Thomas Egger - Mathias Bavay
 * @date   2011-01-02
 * @brief Max range filter.
 * Reject all values greater than the max. Remarks:
 * - the maximum permissible value has to be provided has an argument (in SI)
 * - the keyword "soft" maybe added, in such a case all data greater than the max would be assigned
 * either the maximum permissible value or another value given as an extra argument (350 in the example below)
 * @code
 * TA::filter1	= max
 * TA::arg1	= soft 330 350
 * @endcode
 */

class FilterMax : public FilterBlock {
	public:
		FilterMax(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);

		double max_val;
		double max_soft;
		bool is_soft;
};

} //end namespace

#endif
