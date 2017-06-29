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
#ifndef WINSTRAL_LISTON_ALGORITHM_H
#define WINSTRAL_LISTON_ALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

class WinstralListonAlgorithm : public InterpolationAlgorithm {
	public:
		WinstralListonAlgorithm(Meteo2DInterpolator& i_mi,
		                  const std::vector<std::string>& i_vecArgs,
		                  const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		void initGrid(const DEMObject& dem, Grid2DObject& grid);
		static bool windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station);
		static void getSynopticWind(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station, double& VW, double& DW);

		std::string base_algo, ref_station;
		bool inputIsAllZeroes;
		static const double dmax;
};

} //end namespace mio

#endif
