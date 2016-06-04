/***********************************************************************************/
/*  Copyright 2009 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#ifndef __PGMIO_H__
#define __PGMIO_H__

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/dataClasses/Grid3DObject.h>

#include <string>
#include <sstream>
#include <iostream>

namespace mio {

/**
 * @class PGMIO
 * @brief This class writes 2D grids as Portable Grey Map file format.
 *
 * @ingroup plugins
 * @author Mathias Bavay
 * @date   2010-06-17
 */
class PGMIO : public IOInterface {
	public:
		PGMIO(const std::string& configfile);
		PGMIO(const PGMIO&);
		PGMIO(const Config& cfgreader);
		~PGMIO() throw();

		virtual void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
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

		void read3DGrid(Grid3DObject& grid_out, const std::string& in_name); //HACK

	private:
		void getGridPaths();
		void cleanup() throw();
		void read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name);
		size_t getNextHeader(std::vector<std::string>& vecString, const std::string& filename);

		const Config cfg;
		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::ifstream fin; //Input file streams
		std::ofstream fout;//Output file streams
		std::string grid2dpath_in, grid2dpath_out;
};

} //namespace
#endif
