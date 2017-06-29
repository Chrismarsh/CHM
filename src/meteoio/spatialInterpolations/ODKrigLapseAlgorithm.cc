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

#include <meteoio/spatialInterpolations/ODKrigLapseAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

void LapseOrdinaryKrigingAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	//optimization: getRange (from variogram fit -> exclude stations that are at distances > range (-> smaller matrix)
	//or, get max range from io.ini, build variogram from this user defined max range
	std::vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Fit1D trend(Fit1D::NOISY_LINEAR, vecAltitudes, vecData, false);
	if (!trend.fit())
		throw IOException("Interpolation FAILED for parameter " + MeteoData::getParameterName(param) + ": " + trend.getInfo(), AT);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecData);

	if (!computeVariogram(true)) //only refresh once a month, or once a week, etc
		throw IOException("The variogram for parameter " + MeteoData::getParameterName(param) + " could not be computed!", AT);
	Interpol2D::ODKriging(vecData, vecMeta, dem, variogram, grid);

	retrend(dem, trend, grid);
}

} //namespace
