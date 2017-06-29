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

#include <meteoio/spatialInterpolations/UserAlgorithm.h>
#include <meteoio/FileUtils.h>

namespace mio {

std::string USERInterpolation::getGridFileName() const
{
	const size_t nrArgs = vecArgs.size();
	if (nrArgs > 2) {
		throw InvalidArgumentException("Too many arguments for the "+algo+" interpolation algorithm", AT);
	}
	const std::string prefix = (nrArgs==1)? vecArgs[0] + "/" : "";
	const std::string ext = (nrArgs==2)? vecArgs[1] : ".asc";
	const std::string gridname(  prefix + date.toString(Date::NUM) + "_" + MeteoData::getParameterName(param) + ext );

	return gridname;
}

double USERInterpolation::getQualityRating(const Date& i_date, const MeteoData::Parameters& in_param)
{
	date = i_date;
	param = in_param;
	filename = getGridFileName();

	if (grid2d_path.empty())
		gridsmanager.getConfig().getValue("GRID2DPATH", "Input", grid2d_path);

	if (!FileUtils::validFileAndPath(grid2d_path+"/"+filename)) {
		std::cerr << "[E] Invalid grid filename for "+algo+" interpolation algorithm: " << grid2d_path+"/"+filename << "\n";
		return 0.0;
	}

	return (FileUtils::fileExists(grid2d_path+"/"+filename))? 1. : 0.;
}

void USERInterpolation::calculate(const DEMObject& dem, Grid2DObject& grid)
{
	info.clear(); info.str("");
	gridsmanager.read2DGrid(grid, filename);
	if (!grid.isSameGeolocalization(dem)) {
		throw InvalidArgumentException("[E] trying to load a grid(" + filename + ") that does not have the same georeferencing as the DEM!", AT);
	} else {
		info << FileUtils::getFilename(filename);
	}
}

} //namespace
