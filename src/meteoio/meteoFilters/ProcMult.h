/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __PROCMULT_H__
#define __PROCMULT_H__

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcMult
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2012-02-06
 * @brief Multiply values.
 * This multiplies all values by a given factor. Either a fixed value is given as single argument or a period
 * (hourly/daily/monthly) as well as a filename (and absolute or relative path) containing the factors to apply.
 * This file must contain in the first column the indices (months from 1 to 12 or days from 1 to 366 or hours from 0 to 23)
 * and the matching factor in the second column (<a href="http://www.cplusplus.com/reference/cctype/isspace/">whitespace</a> delimited).
 * Comments following the same syntax as in the ini file are accepted, missing indices are treated as 1.
 * @code
 * PSUM::filter1	= mult
 * PSUM::arg1	= 1.3
 *
 * ISWR::filter1	= mult
 * ISWR::arg1		= monthly input/ISWR_corr.dat
 * @endcode
 *
 * Example of correction file (monthly correction, August will receive a correction of 1):
 * @code
 * 01 0.440593
 * 02 0.815111
 * 03 0.475562
 * 04 0.674975
 * 05 0.700086
 * 06 0.886783
 * 07 1.70733
 * 09 1.26533
 * 10 0.577152
 * 11 0.394095
 * 12 0.347335
 * @endcode
 */

class ProcMult : public ProcessingBlock {
	public:
		ProcMult(const std::vector<std::string>& vec_args, const std::string& name, const std::string& i_root_path);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		void parse_args(const std::vector<std::string>& vec_args);

		std::vector<double> vecFactors;
		std::string root_path;
		double factor;
		char type;
};

} //end namespace

#endif
