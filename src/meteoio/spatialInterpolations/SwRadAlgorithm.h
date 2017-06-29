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
#ifndef SWRADINTERPOLATION_H
#define SWRADINTERPOLATION_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>
#include <meteoio/meteoLaws/Sun.h>

namespace mio {

/**
 * @class SWRadInterpolation
 * @brief %Solar radiation interpolation with optional terrain shading.
 * The splitting coefficients and an atmospheric losses factors are computed at each station that provides ISWR and spatially interpolated
 * with an Inverse Distance Weighting scheme. Then the potential radiation is computed at each pixel and scaled appropriately with the
 * atmospheric loss factor for this pixel. When applying topographic shading (default), the local splitting coefficient is used. The global, horizontal
 * short wave radiation is then returned. To turn off the topographic shading, provide the "no_shading" argument.
 *
 * @code
 * ISWR::algorithms = SWRad
 * ISWR::SWRad = no_shading
 * @endcode
 *
 * @note For this method to work, you also need to define spatial interpolations algorithms for TA, RH and P (a basic STD_PRESS algorithm
 * is usually enough)
 * @note This algorithm is quite time consuming (specially the topographic shading) and therefore not appropriate for very large domains.
 */
class SWRadInterpolation : public InterpolationAlgorithm {
	public:
		SWRadInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
  			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), Sun(), vecIdx(), shading(true) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		SunObject Sun;
		std::vector<size_t> vecIdx;
		bool shading; ///<sould we also compute the shading?
		static const double soil_albedo, snow_albedo, snow_thresh;
};

} //end namespace mio

#endif
