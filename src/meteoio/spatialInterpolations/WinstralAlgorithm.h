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
#ifndef WINSTRAL_ALGORITHM_H
#define WINSTRAL_ALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class WinstralAlgorithm
 * @brief DEM-based wind-exposure interpolation algorithm.
 * This is an implementation of the method described in Winstral, Elder, & Davis,
 * <i>"Spatial snow modeling of wind-redistributed snow using terrain-based parameters"</i>, 2002,
 * Journal of Hydrometeorology, <b>3(5)</b>, 524-538.
 * The DEM is used to compute wind exposure factors that are used to alter the precipitation fields.
 * It is usually a good idea to provide a DEM that also contain the accumulated snow height in order
 * to get a progressive softening of the terrain features.
 *
 * This method must therefore first use another algorithm to generate an initial precipitation field,
 * and then modifies this field accordingly. By default, this base method is "idw_lapse" and switches to
 * "avg" if only one station can provide the precipitation at a given time step.
 *
 * Then it requires a synoptic wind direction that can be provided by different means:
 *  - without any extra argument, the stations are located in the DEM and their wind shading (or exposure)
 * is computed. If at least one station is found that is not sheltered from the wind (in every direction), it
 * provides the synoptic wind (in case of multiple stations, the vector average is used). Please note that
 * the stations that are not included in the DEM are considered to be sheltered. If no such station
 * is found, the vector average of all the available stations is used.
 *  - by providing a fixed synoptic wind bearing that is used for all time steps
 *  - by providing the station_id of the station to get the wind direction from. In this case, the base algorithm
 * for generating the initial wind field must be specified in the first position.
 *
 * @remarks Only cells with an air temperature below freezing participate in the redistribution
 * @code
 * PSUM::algorithms    = WINSTRAL
 * PSUM::winstral = idw_lapse 180
 * @endcode
 */
class WinstralAlgorithm : public InterpolationAlgorithm {
	public:
		WinstralAlgorithm(Meteo2DInterpolator& i_mi,
		                  const std::vector<std::string>& i_vecArgs,
		                  const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		void initGrid(const DEMObject& dem, Grid2DObject& grid);
		static bool windIsAvailable(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station);
		static bool isExposed(const DEMObject& dem, Coords location);
		static double getSynopticBearing(const std::vector<MeteoData>& vecMeteo, const std::string& ref_station);
		static double getSynopticBearing(const std::vector<MeteoData>& vecMeteo);
		static double getSynopticBearing(const DEMObject& dem, const std::vector<MeteoData>& vecMeteo);

		std::string base_algo, ref_station;
		double user_synoptic_bearing;
		bool inputIsAllZeroes;
		static const double dmax;
};

} //end namespace mio

#endif
