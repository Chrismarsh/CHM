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
#ifndef LISTONWINDALGORITHM_H
#define LISTONWINDALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class ListonWindAlgorithm
 * @brief Curvature/slope influenced wind interpolation algorithm.
 * This is an implementation of the method described in G. E. Liston and K. Elder,
 * <i>"A meteorological distribution system for high-resolution terrestrial modeling (MicroMet)"</i>, Journal of Hydrometeorology, <b>7.2</b>, 2006.
 * The wind speed and direction are spatially interpolated using IDWLapseAlgorithm. Then, the wind speed and
 * direction fields are altered by wind weighting factors and wind diverting factors (respectively) calculated
 * from the local curvature and slope (as taken from the DEM, see DEMObject). The wind diverting factor is
 * actually the same as in RyanAlgorithm.
 */
class ListonWindAlgorithm : public InterpolationAlgorithm {
	public:
		ListonWindAlgorithm(Meteo2DInterpolator& i_mi,
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
