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
#ifndef IDWLAPSE_ALGORITHM_H
#define IDWLAPSE_ALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class IDWLapseAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * The input data is detrended and the residuals are spatially interpolated using an Inverse Distance
 * Weighting interpolation algorithm (see IDWAlgorithm). Then, each value is reprojected to the real
 * elevation of the relative cell (re-trending). The lapse rate is either calculated from the data
 * (if no extra argument is provided), or given by the user-provided the optional argument <i>"idw_lapse"</i>.
 * If followed by <i>"soft"</i>, then an attempt to calculate the lapse rate from the data is made, any only if
 * unsuccessful or too bad (r^2<0.6), then the user provided lapse rate is used as a fallback.
 * If the optional user given lapse rate is
 * followed by <i>"frac"</i>, then the lapse rate is understood as a fractional lapse rate, that is a relative change
 * of the value as a function of the elevation (for example, +0.05% per meters given as 0.0005). In this case, no attempt to calculate
 * the fractional lapse from the data is made.
 */
class IDWLapseAlgorithm : public InterpolationAlgorithm {
	public:
		IDWLapseAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		bool lapse_rate_provided; ///< when giving a lapse rate, the requirements on the number of stations are relaxed
};

} //end namespace mio

#endif
