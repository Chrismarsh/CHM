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

#include <meteoio/dataGenerators/SinGenerator.h>
#include <meteoio/MathOptim.h>

namespace mio {

void SinGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if (vecArgs.size()==4) {
		const std::string type_str=IOUtils::strToUpper(vecArgs[0]);
		if ( type_str=="YEARLY" ) type='y';
		else if ( type_str=="DAILY" ) type='d';
		else
			throw InvalidArgumentException("Invalid period \""+type_str+"\" specified for the "+algo+" generator", AT);

		double min, max;
		const bool status1 = IOUtils::convertString(min, vecArgs[1]);
		const bool status2 = IOUtils::convertString(max, vecArgs[2]);
		amplitude = 0.5*(max-min); //the user provides min, max
		offset = min+amplitude;
		const bool status3 = IOUtils::convertString(phase, vecArgs[3]);
		if (!status1 || !status2 || !status3)
			throw InvalidArgumentException("Can not parse the arguments supplied for the "+algo+" generator", AT);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool SinGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value == IOUtils::nodata) {
		double t; //also, the minimum must occur at 0 if phase=0
		if (type=='y') {
			t = (static_cast<double>(md.date.getJulianDayNumber()) - phase*365.25) / 366.25 - .25;
		} else if (type=='d') {
			const double julian = md.date.getJulian();
			t = (julian - Optim::intPart(julian) - phase) + .25; //watch out: julian day starts at noon!
		} else {
			std::ostringstream ss;
			ss << "Invalid period \"" << type << "\" specified for the " << algo << " generator";
			throw InvalidArgumentException(ss.str(), AT);
		}

		static const double w = 2.*Cst::PI;
		value = amplitude * sin(w*t) + offset;
	}

	return true; //all missing values could be filled
}

bool SinGenerator::create(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;

	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		generate(param, vecMeteo[ii]);
	}

	return true; //all missing values could be filled
}

} //namespace
