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
#ifndef __PROCADD_H__
#define __PROCADD_H__

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcAdd
 * @ingroup processing
 * @author Mathias Bavay
 * @date   2012-02-06
 * @brief Add an offset to the values.
 * This adds to all values a given offset. Either a fixed value is given as single argument or a period
 * (hourly/daily/monthly) as well as a filename (and absolute or relative path) containing the offsets to apply.
 * This file must contain in the first column the indices (months from 1 to 12 or days from 1 to 366 or hours from 0 to 23)
 * and the matching offset in the second column (<a href="http://www.cplusplus.com/reference/cctype/isspace/">whitespace</a> delimited).
 * Comments following the same syntax as in the ini file are accepted, missing indices are treated as 0.
 * @code
 * TA::filter1	= add
 * TA::arg1	= 2.5
 *
 * TSG::filter1	= add
 * TSG::arg1	= daily input/TSG_corr.dat
 * @endcode
 *
 * Example of correction file (monthly correction, December will receive a correction of 0):
 * @code
 * 01 -0.375
 * 02 -1.932
 * 03 -4.304
 * 04 -2.449
 * 05 -1.629
 * 06 -1.734
 * 07 -2.414
 * 09 -1.289
 * 10 -1.086
 * 11 -0.769
 * @endcode
 */

class ProcAdd : public ProcessingBlock {
	public:
		ProcAdd(const std::vector<std::string>& vec_args, const std::string& name, const std::string& i_root_path);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	protected:
		static void readCorrections(const std::string& filter, const std::string& filename, const char& c_type, std::vector<double> &corrections);

	private:
		void parse_args(const std::vector<std::string>& vec_args);

		std::vector<double> vecOffsets;
		std::string root_path;
		double offset;
		char type;
};

} //end namespace

#endif
