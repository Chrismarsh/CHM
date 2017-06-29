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
#ifndef AVGLAPSERATEALGORITHM_H
#define AVGLAPSERATEALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class AvgLapseRateAlgorithm
 * @brief Average filling with elevation lapse rate interpolation algorithm.
 * The grid is filled with the average of the detrended measured values and then re-trended. Or to put it
 * differently, the following operations are performed: detrending - averaging - re-trending.
 * The lapse rate is either calculated from the data
 * (if no extra argument is provided), or given by the user-provided the optional argument <i>"avg_lapse"</i>.
 * If followed by <i>"soft"</i>, then an attempt to calculate the lapse rate from the data is made, any only if
 * unsuccessful, then user provided lapse rate is used as a fallback. If the optional user given lapse rate is
 * followed by <i>"frac"</i>, then the lapse rate is understood as a fractional lapse rate, that is a relative change
 * of the value as a function of the elevation (for example, +0.05% per meters given as 0.0005). In this case, no attempt to calculate
 * the fractional lapse from the data is made.
 * @code
 * PSUM::algorithms = AVG_LAPSE
 * PSUM::avg_lapse   = soft 0.05 frac
 * @endcode
 */
class AvgLapseRateAlgorithm : public InterpolationAlgorithm {
	public:
		AvgLapseRateAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
};

} //end namespace mio

#endif
