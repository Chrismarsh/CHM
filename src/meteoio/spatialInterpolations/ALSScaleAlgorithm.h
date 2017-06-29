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
#ifndef ALS_INTERPOLATION_H
#define ALS_INTERPOLATION_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class ALS_Interpolation
 * @brief Scale and distribute the precipitation according to Airborn Laser Scans (ALS) grids.
 * This needs two arguments: first the base method to fill the grid (for example, idw_lapse) and
 * then the name of the file (in GRID2DPATH) containing the gridded ALS data (relying on the GRID2D plugin).
 * If there are some time steps when only one station provides the necessary parameter, the base method will
 * automatically switch to "AVG". A third (optional) argument can be provided that is the air temperature
 * threshold (in K) below which such redistribution occurs (so liquid precipitation is not redistributed).
 *
 * @code
 * PSUM::algorithms = ALS_SCALING
 * PSUM::als_scaling = idw_lapse als_20150213.asc
 * @endcode
 */
class ALS_Interpolation : public InterpolationAlgorithm {
	public:
		ALS_Interpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		void initGrid(const DEMObject& dem, Grid2DObject& grid);
		Grid2DObject ALS_scan;
		std::string filename, grid2d_path, base_algo, base_algo_user;
		double ta_thresh, als_mean; ///< the air temperature must be below a given threshold for the scaling to be applied
		bool inputIsAllZeroes;
};

} //end namespace mio

#endif
