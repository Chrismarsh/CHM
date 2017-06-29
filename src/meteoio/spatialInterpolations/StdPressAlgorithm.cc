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

#include <meteoio/spatialInterpolations/StdPressAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>
#include <meteoio/meteoLaws/Atmosphere.h>

namespace mio {

double StandardPressureAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	const size_t nr_args = vecArgs.size();
	if (nr_args>1)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	if (nr_args==1) {
		if (IOUtils::strToUpper(vecArgs[0])=="USE_RESIDUALS") {
			use_residuals = true;
		} else
			throw InvalidArgumentException("Unknown argument \""+vecArgs[0]+"\" supplied to the "+algo+" interpolation", AT);
	}

	if (param != MeteoData::P) return 0.0;
	if (nrOfMeasurments <=1 || use_residuals) return 1.0;
	return 0.1;
}

void StandardPressureAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	Interpol2D::stdPressure(dem, grid);

	if (nrOfMeasurments==1 || !use_residuals) { //correct the locally measured offset to std pressure
		double offset = 0.;
		size_t count = 0;
		for (size_t ii=0; ii<nrOfMeasurments; ii++) {
			const double altitude = vecMeta[ii].position.getAltitude();
			if (altitude!=IOUtils::nodata) {
				offset += vecData[ii] - Atmosphere::stdAirPressure( altitude );
				count++;
			}
		}
		if (count>0) {
			offset /= static_cast<double>( count );
			grid += offset;
		}
	} else if (use_residuals && nrOfMeasurments>1) { //spatially distribute the residuals
		std::vector<double> residuals;
		for (size_t ii=0; ii<nrOfMeasurments; ii++) {
			const double altitude = vecMeta[ii].position.getAltitude();
			if (altitude!=IOUtils::nodata)
				residuals.push_back( vecData[ii] - Atmosphere::stdAirPressure( altitude ) );
		}
		if (residuals.empty())
			throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

		Grid2DObject offset;
		Interpol2D::IDW(residuals, vecMeta, dem, offset);
		grid += offset;
	}
}

} //namespace
