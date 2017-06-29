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
#ifndef RYANALGORITHM_H
#define RYANALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class RyanAlgorithm
 * @brief DEM-based wind direction interpolation algorithm.
 * This is an implementation of the method described in Ryan,
 * <i>"a mathematical model for diagnosis and prediction of surface winds in mountainous terrain"</i>,
 * 1977, journal of applied meteorology, <b>16</b>, 6.
 * The DEM is used to compute wind drection changes that are used to alter the wind direction fields.
 * @code
 * DW::algorithms    = RYAN
 * @endcode
 */
class RyanAlgorithm : public InterpolationAlgorithm {
	public:
		RyanAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), vecDataVW(), vecDataDW(), inputIsAllZeroes(false) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataVW, vecDataDW; ///<vectors of extracted VW and DW
		bool inputIsAllZeroes;
};

} //end namespace mio

#endif
