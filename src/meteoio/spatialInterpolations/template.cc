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

#include <meteoio/spatialInterpolations/template.h>

namespace mio {

double TEMPLATE::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	nrOfMeasurments = getData(date, param, vecData, vecMeta);

	//parse the arguments, for example:
	const size_t nr_args = vecArgs.size();
	if (nr_args>1)
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" algorithm", AT);
	if (nr_args==1) {
		//parse one argument
	}

	//depending on the arguments and the available data, return a quality rating.
	//Have a look at the other method in order to figure out which rating should be return to rank
	//where you want within the other algorithms
	return 0.1;
}

void TEMPLATE::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	//this should contain all the information to forward to the user so he can understand how the interpolation went
	info.clear(); info.str("");

	//depending on how the interpolation is done, either reset the grid here or in a call to a method of  Interpol2D
	grid.set(dem, IOUtils::nodata);

	//now fill each cell of the grid

}

} //namespace
