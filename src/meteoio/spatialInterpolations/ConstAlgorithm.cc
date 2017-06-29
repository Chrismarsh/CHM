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

#include <meteoio/spatialInterpolations/ConstAlgorithm.h>
#include <meteoio/meteoStats/libinterpol2D.h>

namespace mio {

double ConstAlgorithm::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	const size_t nr_args = vecArgs.size();
	if (nr_args!=1)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);

	IOUtils::convertString(user_cst, vecArgs[0]);
	return 0.01;
}

void ConstAlgorithm::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	Interpol2D::constant(user_cst, dem, grid);
}

} //namespace
