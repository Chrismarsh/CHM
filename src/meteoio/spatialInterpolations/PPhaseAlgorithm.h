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
#ifndef PPHASEINTERPOLATION_H
#define PPHASEINTERPOLATION_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class PPHASEInterpolation
 * @brief Precipitation phase splitting generation
 * This does not interpolate any measured precipitation phase but generates it for each point based on parametrizations, similarly to the PPHASE generator
 * (see PPhaseGenerator).
 *
 * The methods that are offered are currently the following:
 * - THRESH: a provided fixed air temperature threshold splits precipitation as either fully solid or fully liquid
 * - RANGE: two air temperature thresholds provide the lower and upper range for fully solid / fully liquid precipitation.
 *                 Within the provided range, a linear transition is assumed.
 * @code
 * PSUM::algorithms = PPHASE
 * PSUM::pphase = THRESH 274.35
 * @endcode
 */
class PPHASEInterpolation : public InterpolationAlgorithm {
	public:
		PPHASEInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
  			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager),
  			model(THRESH), fixed_thresh(IOUtils::nodata), range_start(IOUtils::nodata), range_norm(IOUtils::nodata) {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		typedef enum PARAMETRIZATION {
				THRESH,
				RANGE
			} parametrization;
		parametrization model;
		double fixed_thresh, range_start, range_norm;
};

} //end namespace mio

#endif
