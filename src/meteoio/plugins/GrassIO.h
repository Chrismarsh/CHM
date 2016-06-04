/***********************************************************************************/
/*  Copyright 2009 EPFL                                                            */
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
#ifndef __GRASSIO_H__
#define __GRASSIO_H__

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace mio {

/**
 * @class GrassIO
 * @brief This class enables the access to 2D grids stored in GRASS ASCII (e.g. JGrass) format
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2008-08-03
 */
class GrassIO : public IOInterface {
	public:
		GrassIO(const std::string& configfile);
		GrassIO(const GrassIO&);
		GrassIO(const Config&);
		~GrassIO() throw();

		virtual void read2DGrid(Grid2DObject& dem_out, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);

		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);

		virtual void writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo,
		                            const std::string& name="");

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readPOI(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& filename);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void cleanup() throw();

		const Config cfg;
		std::ifstream fin; //Input file streams
		std::ofstream fout;//Output file streams
		static const double plugin_nodata;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
};

} //end namespace mio

#endif
