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
#ifndef LOCALIDWLAPSE_ALGORITHM_H
#define LOCALIDWLAPSE_ALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class LocalIDWLapseAlgorithm
 * @brief Inverse Distance Weighting interpolation algorithm with elevation detrending/reprojection.
 * The closest n stations (n being given as an extra argument of <i>"lidw_lapse"</i>) to each pixel are
 * used to compute the local lapse rate, allowing to project the contributions of these n stations to the
 * local pixel with an inverse distance weight. Beware, this method sometimes produces very sharp transitions
 * as it spatially moves from one station's area of influence to another one!
 */
class LocalIDWLapseAlgorithm : public InterpolationAlgorithm {
	public:
		LocalIDWLapseAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		size_t nrOfNeighbors;
};

} //end namespace mio

#endif
