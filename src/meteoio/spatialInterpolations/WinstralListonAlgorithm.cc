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

#include <meteoio/spatialInterpolations/WinstralListonAlgorithm.h>
#include <meteoio/spatialInterpolations/WinstralAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

const double WinstralListonAlgorithm::dmax = 300.;

WinstralListonAlgorithm::WinstralListonAlgorithm(Meteo2DInterpolator& i_mi, const std::vector<std::string>& i_vecArgs,
                                     const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
                  : InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), base_algo("IDW_LAPSE"), ref_station(),
                    inputIsAllZeroes(false)
{
	const size_t nr_args = vecArgs.size();
	if (nr_args==2) {
		base_algo = IOUtils::strToUpper( vecArgs[0] );
		ref_station = vecArgs[1];
		return;
	} else if (nr_args!=2)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
}

double WinstralListonAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	//This algorithm is only valid for PSUM (we could add HS later)
	if (in_param!=MeteoData::PSUM)
		return 0.0;

	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);
	inputIsAllZeroes = Interpol2D::allZeroes(vecData);

	if (inputIsAllZeroes) return 0.99;

	if (nrOfMeasurments==0) return 0.0;

	if (nrOfMeasurments==1 && ref_station.empty()) { //ie: still using default base_algo
		base_algo = "AVG";
	}

	//check that the necessary wind data is available
	if (!ref_station.empty()) {
		if (!windIsAvailable(vecMeteo, ref_station))
			return 0.0;
	}

	return 0.99;
}

void WinstralListonAlgorithm::initGrid(const DEMObject& dem, Grid2DObject& grid)
{
	//initialize precipitation grid with user supplied algorithm (IDW_LAPSE by default)
	std::vector<std::string> vecArgs2;
	mi.getArgumentsForAlgorithm(MeteoData::getParameterName(param), base_algo, vecArgs2);
	std::auto_ptr<InterpolationAlgorithm> algorithm(AlgorithmFactory::getAlgorithm(base_algo, mi, vecArgs2, tsmanager, gridsmanager));
	algorithm->getQualityRating(date, param);
	algorithm->calculate(dem, grid);
	info << algorithm->getInfo();
}

bool WinstralListonAlgorithm::windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station)
{
	if (ref_station.empty()) {
		for (size_t ii=0; ii<vecMeteo.size(); ii++) {
			const double VW = vecMeteo[ii](MeteoData::VW);
			const double DW = vecMeteo[ii](MeteoData::DW);
			if (VW!=IOUtils::nodata && DW!=IOUtils::nodata)
				return true; //at least one station is enough
		}
	} else {
		double VW, DW;
		getSynopticWind(vecMeteo, ref_station, VW, DW);
		if (VW!=IOUtils::nodata && DW!=IOUtils::nodata)
			return true;
	}

	return false;
}

void WinstralListonAlgorithm::getSynopticWind(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station, double& VW, double& DW)
{
	for (size_t ii=0; ii<vecMeteo.size(); ++ii) {
		if (vecMeteo[ii].meta.stationID==ref_station) {
			VW = vecMeteo[ii](MeteoData::VW);
			DW = vecMeteo[ii](MeteoData::DW);
			return;
		}
	}
}

void WinstralListonAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");

	//if all data points are zero, simply fill the grid with zeroes
	if (inputIsAllZeroes) {
		Interpol2D::constant(0., dem, grid);
		return;
	}

	//get initial PSUM grid
	initGrid(dem, grid);

	//get meteo fields interpolation from call back to Meteo2DInterpolator
	Grid2DObject ta, dw, vw;
	mi.interpolate(date, dem, MeteoData::TA, ta);
	mi.interpolate(date, dem, MeteoData::DW, dw);
	mi.interpolate(date, dem, MeteoData::VW, vw);

	//alter the field with Winstral and the chosen wind direction
	Interpol2D::Winstral(dem, ta, dw, vw, dmax, grid);
}

} //namespace
