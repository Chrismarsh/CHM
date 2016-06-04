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
#ifndef __BORMAIO_H__
#define __BORMAIO_H__

#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>
#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/Date.h>

#include <libxml++/libxml++.h>
#include <string>
#include <sstream>
#include <iostream>

namespace mio {

/**
 * @class BormaIO
 * @brief This class enables the access meteo data in Borma's XML format
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2008-11-20
 */

class BormaIO : public IOInterface {
	public:
		BormaIO(void (*delObj)(void*), const Config& i_cfg);

		BormaIO(const std::string& configfile);
		BormaIO(const BormaIO&);
		BormaIO(const Config&);

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

		virtual void write2DGrid(const Grid2DObject& grid_in, const std::string& name);
		virtual void write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date);

	private:
		void convertUnits(MeteoData& meteo);
		void checkForMeteoFiles(const std::string& xmlpath, const std::string& stationname, const Date& date_in,
		                        std::string& filename_out, Date& date_out);
		void xmlParseStringToDouble(const std::string& str_in, double& d_out, const std::string& parname);
		std::string xmlGetNodeContent(xmlpp::Node* pNode, const std::string& nodename);
		void xmlExtractData(const std::string& filename, const Date& date_in, MeteoData& md, StationData& sd);
		std::string xmlGetNodeName(xmlpp::Node* pNode);
		xmlpp::Node* xmlGetNode(xmlpp::Node* parentNode, const std::string& nodename);
		Date stringToDate(const std::string& tmp) const;
		bool validFilename(const std::string& tmp) const;
		void getFiles(const std::string& stationsname, const Date& start_date, const Date& end_date,
		              std::vector<std::string>& vecFiles, std::vector<Date>& vecDate);
		void readStationNames(void);
		bool bufferData(const Date& dateStart, const Date& dateEnd,
		                std::vector< std::vector<MeteoData> >& vecMeteo,
		                const unsigned int& stationnr);

		std::vector<std::string> vecStationName;
		const Config cfg;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		double in_tz;
		size_t nr_stations; //number of stations to read from

		static const double plugin_nodata; //plugin specific nodata value, e.g. -999
		static const double default_tz; //default timezone
		static const double pivot_year; //pivot year for Y2K suppport
		static const std::string dflt_extension;


};

} //end namespace mio

#endif
