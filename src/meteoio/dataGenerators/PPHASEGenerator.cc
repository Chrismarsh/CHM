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

#include <meteoio/dataGenerators/PPHASEGenerator.h>

namespace mio {

void PPhaseGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	const size_t nArgs = vecArgs.size();
	
	if (nArgs<1 || IOUtils::isNumeric(vecArgs[0]))
		throw InvalidArgumentException("Wrong arguments supplied to the "+algo+" generator. Please provide the method to use and its arguments!", AT);
	
	const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);
	if (user_algo=="THRESH") {
		if (nArgs!=2)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" generator for the "+user_algo+" method", AT);
		const bool status = IOUtils::convertString(fixed_thresh, vecArgs[1]);
		if (!status)
			throw InvalidArgumentException(algo+" generator, "+user_algo+" method: can not parse provided threshold", AT);
		model = THRESH;
	} else if (user_algo=="RANGE") {
		if (nArgs!=3)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" generator for the "+user_algo+" method", AT);
		double range_thresh1, range_thresh2;
		const bool status1 = IOUtils::convertString(range_thresh1, vecArgs[1]);
		const bool status2 = IOUtils::convertString(range_thresh2, vecArgs[2]);
		if (!status1 || !status2)
			throw InvalidArgumentException(algo+" generator, "+user_algo+" method: can not parse provided thresholds", AT);
		if (range_thresh1==range_thresh2)
			throw InvalidArgumentException(algo+" generator, "+user_algo+" method: the two provided threshold must be different", AT);
		if (range_thresh1>range_thresh2) 
			std::swap(range_thresh1, range_thresh2);
		range_start = range_thresh1;
		range_norm = 1. / (range_thresh2-range_thresh1);
		model = RANGE;
	} else
		throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied to the "+algo+" generator", AT);
}

bool PPhaseGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA);
		if (TA==IOUtils::nodata) return false;
		
		if (model==THRESH) {
			value = (TA>=fixed_thresh)? 1. : 0.;
		} else if (model==RANGE) {
			const double tmp_rainfraction = range_norm * (TA - range_start);
			value = (tmp_rainfraction>1)? 1. : (tmp_rainfraction<0.)? 0. : tmp_rainfraction;
		}
	}

	return true; //all missing values could be filled
}

bool PPhaseGenerator::create(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;
	
	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (!generate(param, vecMeteo[ii]))
			all_filled = false;
	}

	return all_filled;
}

} //namespace
