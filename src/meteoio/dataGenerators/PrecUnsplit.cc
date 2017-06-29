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

#include <meteoio/dataGenerators/PrecUnsplit.h>

namespace mio {

bool PrecUnsplit::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value == IOUtils::nodata) {
		if (!md.param_exists("RAINF") || !md.param_exists("SNOWF")) return false;

		const double rain = md("RAINF");
		const double snow = md("SNOWF");
		if (rain==IOUtils::nodata || snow==IOUtils::nodata) return false;

		const double PSUM = rain+snow;

		if (param==MeteoData::PSUM) {
			md(MeteoData::PSUM) = PSUM;
		} else {
			if (PSUM==0.)
				md(MeteoData::PSUM_PH) = 1.;
			else
				md(MeteoData::PSUM_PH) = rain / PSUM;
		}
	}

	return true; //all missing values could be filled
}

bool PrecUnsplit::create(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		const bool status = generate(param, vecMeteo[ii]);
		if (status==false) all_filled = false;
	}

	return all_filled;
}

} //namespace
