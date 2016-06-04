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
#ifndef __IOHANDLER_H__
#define __IOHANDLER_H__

#include <meteoio/IOInterface.h>
#include <meteoio/IOExceptions.h>

#include <map>
#include <set>
#include <string>

namespace mio {

/**
* @file IOHandler.h
* @class IOHandler
* @brief This class is the class to use for raw I/O operations. It is responsible for transparently loading the plugins
* and it follows the interface defined by the IOInterface class with the addition of a few convenience methods.
*/
class IOHandler : public IOInterface {
	public:
		IOHandler(const IOHandler&);
		IOHandler(const Config&);

		virtual ~IOHandler() throw();

		IOHandler& operator=(const IOHandler&); ///<Assignement operator

		//methods defined in the IOInterface class
		virtual void read2DGrid(Grid2DObject& out_grid, const std::string& parameter="");
		virtual void read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date);
		virtual void readDEM(DEMObject& dem_out);
		virtual void readLanduse(Grid2DObject& landuse_out);
		virtual void readStationData(const Date& date,
		                             STATIONS_SET& vecStation);

		virtual void writeMeteoData(const std::vector<METEO_SET>& vecMeteo,
		                            const std::string& name="");
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector<METEO_SET>& vecMeteo,
		                           const size_t& stationindex=IOUtils::npos);

		virtual void readAssimilationData(const Date&, Grid2DObject& da_out);
		virtual void readPOI(std::vector<Coords>& pts);
		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

		const std::string toString() const;

	private:
		IOInterface* getPlugin(const std::string& plugin_name) const;
		IOInterface* getPlugin(const std::string& cfgkey, const std::string& cfgsection);
		void parse_copy_config();
		void create_exclude_map();
		void create_keep_map();
		void checkTimestamps(const std::vector<METEO_SET>& vecVecMeteo) const;
		void exclude_params(std::vector<METEO_SET>& vecVecMeteo) const;
		void keep_params(std::vector<METEO_SET>& vecVecMeteo) const;
		void copy_parameters(const size_t& stationindex, std::vector< METEO_SET >& vecMeteo) const;

		const Config& cfg;
		std::map<std::string, IOInterface*> mapPlugins;
		std::map< std::string, std::set<std::string> > excluded_params;
		std::map< std::string, std::set<std::string> > kept_params;
		std::vector<std::string> copy_parameter, copy_name;
		bool enable_copying, excludes_ready, keeps_ready;
};

} //namespace

#endif
