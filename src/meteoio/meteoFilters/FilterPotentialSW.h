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
#ifndef FILTERPOTENTIALSW_H
#define FILTERPOTENTIALSW_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  FilterPotentialSW
 * @ingroup processing
 * @brief Checks for physically reallistic incoming short wave radiation (ISWR) values.
 * For each data point, the measured value must be:
 *     + more than the horizontal top of atmosphere potential radiation multiplied by the *min_coeff* coefficient;
 *     + less than the global horizontal potential radiation multiplied by the *max_coeff* coefficient.
 *
 * Either no arguments are provided and these coefficients have the default values of 0.03 and 1.1 respectively or both
 * must be provided. The default values come from Moradi, I., <i>"Quality control of global solar radiation using
 * sunshine duration hours"</i>, 2009, Energy 34, <b>no. 1</b>, 1-6.
 *
 * @code
 * ISWR::filter1	= PotentialSW
 * ISWR::arg1	= 0.03 1.1
 * @endcode
 */

class FilterPotentialSW : public FilterBlock {
	public:
		FilterPotentialSW(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(std::vector<std::string> vec_args);

		double min_coeff, max_coeff;
};

} //end namespace

#endif
