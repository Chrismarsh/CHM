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
#ifndef __ARCIO_H__
#define __ARCIO_H__

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>
#include <iostream>

namespace mio {

/**
 * @class ARCIO
 * @brief This class enables the access to 2D grids stored in ESRI ASCII (ARCGIS) format
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2009-08-28
 */

class ARCIO : public IOInterface {
	public:
		ARCIO(const std::string& configfile);
		ARCIO(const ARCIO&);
		ARCIO(const Config&);
		~ARCIO() throw();

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
		void getGridPaths();
		void cleanup() throw();
		void read2DGrid_internal(Grid2DObject& grid_out, const std::string& full_name);
		const Config cfg;

		std::ifstream fin; //Input file streams
		std::ofstream fout;//Output file streams
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::string grid2dpath_in, grid2dpath_out;
		std::string grid2d_ext_in, grid2d_ext_out; //file extension

		bool a3d_view_in, a3d_view_out; ///< make filename compatible with the Alpine3D's viewer?
};

} //end namespace mio

#endif
