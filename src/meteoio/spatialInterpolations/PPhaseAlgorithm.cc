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

#include <meteoio/spatialInterpolations/PPhaseAlgorithm.h>

namespace mio {

double PPHASEInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;

	const size_t nArgs = vecArgs.size();

	if (nArgs<1 || IOUtils::isNumeric(vecArgs[0]))
		throw InvalidArgumentException("Wrong arguments supplied to the "+algo+" interpolation. Please provide the method to use and its arguments!", AT);

	const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);
	if (user_algo=="THRESH") {
		if (nArgs!=2)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" interpolation for the "+user_algo+" method", AT);
		IOUtils::convertString(fixed_thresh, vecArgs[1]);
		model = THRESH;
	} else if (user_algo=="RANGE") {
		if (nArgs!=3)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" interpolation for the "+user_algo+" method", AT);
		double range_thresh1, range_thresh2;
		IOUtils::convertString(range_thresh1, vecArgs[1]);
		IOUtils::convertString(range_thresh2, vecArgs[2]);
		if (range_thresh1==range_thresh2)
			throw InvalidArgumentException(algo+" interpolation, "+user_algo+" method: the two provided threshold must be different", AT);
		if (range_thresh1>range_thresh2)
			std::swap(range_thresh1, range_thresh2);
		range_start = range_thresh1;
		range_norm = 1. / (range_thresh2-range_thresh1);
		model = RANGE;
	} else
		throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied to the "+algo+" interpolation", AT);

	return 0.1;
}

void PPHASEInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator

	grid.set(dem, IOUtils::nodata);

	if (model==THRESH) {
		for (size_t ii=0; ii<dem.size(); ii++) {
			const double TA=ta(ii);
			if (TA==IOUtils::nodata) continue;
			grid(ii) = (TA>=fixed_thresh)? 1. : 0.;
		}
	} else if (model==RANGE) {
		for (size_t ii=0; ii<dem.size(); ii++) {
			const double TA=ta(ii);
			if (TA==IOUtils::nodata) continue;
			const double tmp_rainfraction = range_norm * (TA - range_start);
			grid(ii) = (tmp_rainfraction>1)? 1. : (tmp_rainfraction<0.)? 0. : tmp_rainfraction;
		}
	}
}

} //namespace
