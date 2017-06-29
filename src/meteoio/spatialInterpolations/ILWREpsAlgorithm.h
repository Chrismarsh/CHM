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
#ifndef ILWREPSALGORITHM_H
#define ILWREPSALGORITHM_H

#include <meteoio/spatialInterpolations/InterpolationAlgorithms.h>

namespace mio {

/**
 * @class ILWREpsAlgorithm
 * @brief Incoming Long Wave Radiation interpolation algorithm.
 * Each ILWR is converted to an emissivity (using the local air temperature), interpolated using AVG_LAPSE or IDW_LAPSE with
 * a fixed lapse rate and reconverted to ILWR.
 *
 * As a side effect, the user must have defined algorithms to be used for air temperature (since this is needed for
 * emissivity to ILWR conversion)
 */
class ILWREpsAlgorithm : public InterpolationAlgorithm {
	public:
		ILWREpsAlgorithm(Meteo2DInterpolator& i_mi,
					const std::vector<std::string>& i_vecArgs,
					const std::string& i_algo, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
			: InterpolationAlgorithm(i_mi, i_vecArgs, i_algo, i_tsmanager, i_gridsmanager), vecDataEA() {}
		virtual double getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param);
		virtual void calculate(const DEMObject& dem, Grid2DObject& grid);
	private:
		std::vector<double> vecDataEA; ///<vectors of extracted emissivities
};

} //end namespace mio

#endif
