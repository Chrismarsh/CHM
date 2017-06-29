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
#ifndef GRIDSMANAGER_H
#define GRIDSMANAGER_H

#include <meteoio/dataClasses/Buffer.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/IOHandler.h>
#include <meteoio/Config.h>

namespace mio {

class GridsManager {
	public:
		GridsManager(IOHandler& in_iohandler, const Config& in_cfg);

		//Legacy support to support functionality of the IOInterface superclass:
		void read2DGrid(Grid2DObject& grid_out, const std::string& parameter="");
		void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		void read3DGrid(Grid3DObject& grid_out, const std::string& i_filename="");
		void read3DGrid(Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		void readDEM(DEMObject& dem_out);
		void readAssimilationData(const Date& date_in, Grid2DObject& da_out);
		void readLanduse(Grid2DObject& landuse_out);
		void write2DGrid(const Grid2DObject& grid_in, const std::string& options="");
		void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);
		void write3DGrid(const Grid3DObject& grid_out, const std::string& options="");
		void write3DGrid(const Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		//end legacy support

		void setProcessingLevel(const unsigned int& i_level);
		void clear_cache();
		
		/**
		 * @brief Returns a copy of the internal Config object.
		 * This is convenient to clone an iomanager
		 * @return new Config object as a copy of the internal Config
		 */
		const Config getConfig() const {return cfg;}

		/**
		 * @brief Returns a copy of the internal IOHandler object.
		 * This is convenient to clone an iomanager
		 * @return new IOHandler object as a copy of the internal IOHandler
		 */
		IOHandler& getIOHandler() const {return iohandler;}

		const std::string toString() const;

	private:
		void addToBuffer(const Grid2DObject& in_grid2Dobj, const std::string& grid_hash);
		bool getFromBuffer(const std::string& grid_hash, Grid2DObject& grid) const;

		IOHandler& iohandler;
		const Config& cfg;
		GridBuffer buffer;

		unsigned int processing_level;
};
} //end namespace
#endif
