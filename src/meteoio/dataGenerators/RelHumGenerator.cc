/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#include <meteoio/dataGenerators/RelHumGenerator.h>
#include <meteoio/meteoLaws/Atmosphere.h>

namespace mio {

bool RhGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value == IOUtils::nodata) {
		const double TA = md(MeteoData::TA);
		if (TA==IOUtils::nodata) return false;//nothing else we can do here

		//first chance to compute RH
		if (md.param_exists("TD")) {
			const double TD = md("TD");
			if (TD!=IOUtils::nodata) {
				value = Atmosphere::DewPointtoRh(TD, TA, false);
				return true;
			}
		}

		//second chance to try to compute RH
		if (md.param_exists("SH")) {
			const double SH = md("SH");
			const double altitude = md.meta.position.getAltitude();
			if (SH!=IOUtils::nodata && altitude!=IOUtils::nodata) {
				value = Atmosphere::specToRelHumidity(altitude, TA, SH);
				return true;
			}
		}

		return false;
	}

	return true; //all missing values could be filled
}

bool RhGenerator::create(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	const double altitude = vecMeteo.front().meta.position.getAltitude(); //if the stations move, this has to be in the loop

	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if (value == IOUtils::nodata) {
			const double TA = vecMeteo[ii](MeteoData::TA);
			if (TA==IOUtils::nodata) { //nothing else we can do here
				all_filled=false;
				continue;
			}

			//first chance to compute RH
			if (vecMeteo[ii].param_exists("TD")) {
				const double TD = vecMeteo[ii]("TD");
				if (TD!=IOUtils::nodata) {
					value = Atmosphere::DewPointtoRh(TD, TA, false);
					continue;
				}
			}

			//second chance to try to compute RH
			if (vecMeteo[ii].param_exists("SH")) {
				const double SH = vecMeteo[ii]("SH");
				if (SH!=IOUtils::nodata && altitude!=IOUtils::nodata) {
					value = Atmosphere::specToRelHumidity(altitude, TA, SH);
					continue;
				}
			}

			all_filled=false;
		}
	}

	return all_filled;
}

} //namespace
