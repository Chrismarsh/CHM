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

#include <meteoio/spatialInterpolations/IDWLapseAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

IDWLapseAlgorithm::IDWLapseAlgorithm(Meteo2DInterpolator& i_mi, const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), lapse_rate_provided(false)
{
	const size_t nr_args = vecArgs.size();
	if (nr_args>2)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	if (nr_args>=1) {
		if ( IOUtils::isNumeric(vecArgs[0])) lapse_rate_provided = true;
	}
}

double IDWLapseAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	if (nrOfMeasurments == 0) return 0.0;
	if (!lapse_rate_provided && nrOfMeasurments<2) return 0.0;

	return 0.7;
}

void IDWLapseAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	std::vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Fit1D trend;
	getTrend(vecAltitudes, vecData, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecData);
	Interpol2D::IDW(vecData, vecMeta, dem, grid); //the meta should NOT be used for elevations!
	retrend(dem, trend, grid);
}

} //namespace
