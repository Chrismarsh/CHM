/***********************************************************************************/
/*  Copyright 2017 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef PrecUnsplit_H
#define PrecUnsplit_H

#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

namespace mio {

/**
 * @class PrecUnsplit
 * @brief This generator converts split precipitation (as provided by the "RAINF" and "SNOWF" parameters)  into amount/phase. It therefore
 * builds either PSUM (if it is attached to PSUM) or PSUM_PH (for any other parameter it is attached to). Please note that
 * when either "RAINF" or "SNOWF" is nodata, no value is generated.
 *
 * @code
 * PSUM_PH::generators = PrecUnsplit
 * PSUM::generators = PrecUnsplit
 * @endcode
 */
class PrecUnsplit : public GeneratorAlgorithm {
	public:
		PrecUnsplit(const std::vector<std::string>& vecArgs, const std::string& i_algo)
			: GeneratorAlgorithm(vecArgs, i_algo) {}
		bool generate(const size_t& param, MeteoData& md);
		bool create(const size_t& param, std::vector<MeteoData>& vecMeteo);
};

} //end namespace mio

#endif
