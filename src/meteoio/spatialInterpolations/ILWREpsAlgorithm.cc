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

#include <meteoio/spatialInterpolations/ILWREpsAlgorithm.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/meteoStats/libinterpol2D.h>
#include <meteoio/meteoLaws/Atmosphere.h>

namespace mio {

double ILWREpsAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for ILWR
	if (in_param != MeteoData::ILWR)
		return 0.0;

	date = i_date;
	param = in_param;
	vecData.clear(); vecMeta.clear();
	vecDataEA.clear();

	nrOfMeasurments = 0;
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		if ((vecMeteo[ii](MeteoData::ILWR) != IOUtils::nodata) && (vecMeteo[ii](MeteoData::TA) != IOUtils::nodata)){
			vecDataEA.push_back( Atmosphere::blkBody_Emissivity( vecMeteo[ii](MeteoData::ILWR), vecMeteo[ii](MeteoData::TA)) );
			vecMeta.push_back(vecMeteo[ii].meta);
			nrOfMeasurments++;
		}
	}

	if (nrOfMeasurments==0)
		return 0.0;

	return 0.9;
}

void ILWREpsAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	std::vector<double> vecAltitudes;
	getStationAltitudes(vecMeta, vecAltitudes);
	if (vecAltitudes.empty())
		throw IOException("Not enough data for spatially interpolating parameter " + MeteoData::getParameterName(param), AT);

	Grid2DObject ta;
	mi.interpolate(date, dem, MeteoData::TA, ta); //get TA interpolation from call back to Meteo2DInterpolator

	Fit1D trend;
	getTrend(vecAltitudes, vecDataEA, trend);
	info << trend.getInfo();
	detrend(trend, vecAltitudes, vecDataEA);
	Interpol2D::IDW(vecDataEA, vecMeta, dem, grid); //the meta should NOT be used for elevations!
	retrend(dem, trend, grid);

	//Recompute Rh from the interpolated td
	for (size_t ii=0; ii<grid.size(); ii++) {
		double &value = grid(ii);
		value = Atmosphere::blkBody_Radiation(value, ta(ii));
	}
}

} //namespace
