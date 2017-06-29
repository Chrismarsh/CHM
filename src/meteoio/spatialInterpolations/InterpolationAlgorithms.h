/***********************************************************************************/
/*  Copyright 2010 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef INTERPOLATIONALGORITHMS_H
#define INTERPOLATIONALGORITHMS_H

#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/TimeSeriesManager.h>
#include <meteoio/GridsManager.h>
#include <meteoio/Meteo2DInterpolator.h>
#include <meteoio/meteoStats/libfit1D.h>

#include <vector>
#include <string>

namespace mio {

class Meteo2DInterpolator; // forward declaration, cyclic header include

/**
 * @class InterpolationAlgorithm
 * @brief A class to perform 2D spatial interpolations. For more, see \ref interpol2d
 *
 * @ingroup stats
 * @author Thomas Egger
 * @date   2010-04-01
*/
class InterpolationAlgorithm {

	public:
		InterpolationAlgorithm(Meteo2DInterpolator& i_mi,
		                       const std::vector<std::string>& i_vecArgs,
		                       const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager) :
		                      algo(i_algo), mi(i_mi), tsmanager(i_tsmanager), gridsmanager(i_gridsmanager), date(0., 0), vecArgs(i_vecArgs), vecMeteo(), vecData(),
		                      vecMeta(), info(), param(MeteoData::firstparam), nrOfMeasurments(0) {}
		virtual ~InterpolationAlgorithm() {}
		//if anything is not ok (wrong parameter for this algo, insufficient data, etc) -> return zero
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param) = 0;
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid) = 0;
		std::string getInfo() const;
		const std::string algo;

 	protected:
		size_t getData(const Date& i_date, const MeteoData::Parameters& i_param, std::vector<double>& o_vecData);
		size_t getData(const Date& i_date, const MeteoData::Parameters& i_param,
		               std::vector<double>& o_vecData, std::vector<StationData>& o_vecMeta);
		static size_t getStationAltitudes(const std::vector<StationData>& i_vecMeta, std::vector<double>& o_vecData);
		void getTrend(const std::vector<double>& vecAltitudes, const std::vector<double>& vecDat, Fit1D &trend) const;
		static void detrend(const Fit1D& trend, const std::vector<double>& vecAltitudes, std::vector<double> &vecDat, const double& min_alt=-1e4, const double& max_alt=1e4);
		static void retrend(const DEMObject& dem, const Fit1D& trend, Grid2DObject &grid, const double& min_alt=-1e4, const double& max_alt=1e4);
		void simpleWindInterpolate(const DEMObject& dem, const std::vector<double>& vecDataVW, const std::vector<double>& vecDataDW, Grid2DObject &VW, Grid2DObject &DW);

		Meteo2DInterpolator& mi;
		TimeSeriesManager& tsmanager;
		GridsManager& gridsmanager;
		Date date;
		const std::vector<std::string> vecArgs; //we must keep our own copy, it is different for each algorithm!

		std::vector<MeteoData> vecMeteo;
		std::vector<double> vecData; ///<store the measurement for the given parameter
		std::vector<StationData> vecMeta; ///<store the station data for the given parameter
		std::ostringstream info; ///<to store some extra information about the interplation process
		MeteoData::Parameters param; ///<the parameter that we will interpolate
		size_t nrOfMeasurments; ///<the available number of measurements
};

class AlgorithmFactory {
	public:
		static InterpolationAlgorithm* getAlgorithm(const std::string& i_algoname,
		                                            Meteo2DInterpolator& i_mi,
		                                            const std::vector<std::string>& i_vecArgs, TimeSeriesManager& tsm, GridsManager& gdm);
};

} //end namespace mio

#endif
