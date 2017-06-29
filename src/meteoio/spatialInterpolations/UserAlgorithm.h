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
#ifndef USERINTERPOLATION_H
#define USERINTERPOLATION_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class USERInterpolation
 * @brief Reads user provided gridded data on the disk.
 * The grids are all in the GRID2DPATH directory given in the [Input] section or in one of
 * its sub-directories that is given as the algorithm's argument (optional). By default, the file extension is assumed to
 * be ".asc" but it is possible to provide as second argument another file extension (then it is mandatory to
 * also provide a sub-directory argument in first position).
 * The files must be named according to the following schema: <b>{numeric date with second resolution}_{capitalized meteo parameter}.{ext}</b>, for example 20081201150000_TA.asc
 * The meteo parameters can be found in \ref meteoparam "MeteoData". Example of use:
 * @code
 * TA::algorithms = USER	# read grids from GRID2DPATH using the GRID2D plugin
 *
 * VW::algorithms = USER	# read grids from GRID2DPATH/wind
 * VW::user       = wind
 *
 * HNW::algorithms = USER	# read grids from GRID2DPATH/precip with the ".dat" extension
 * HNW::user       = precip .dat
 * @endcode
 *
 * If no grid exists for a given timestamp and parameter, the algorithm returns a zero rating so any other interpolation algorithm can pickup
 * and provide a fallback. Therefore, it is not necessary to provide grids for all time steps but one can focuss on only the relevant and interesting
 * time steps.
 */
class USERInterpolation : public InterpolationAlgorithm {
	public:
		USERInterpolation(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), filename(), grid2d_path() {nrOfMeasurments=0;}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::string getGridFileName() const;
		std::string filename, grid2d_path;
};

} //end namespace mio

#endif
