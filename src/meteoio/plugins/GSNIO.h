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
#ifndef __GSNIO_H__
#define __GSNIO_H__


#include <meteoio/Config.h>
#include <meteoio/IOInterface.h>
#include <meteoio/IOUtils.h>

#include <string>
#include <vector>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {

/**
 * @class GSNIO
 * @brief This class enables the access to the GSN RESTful web service
 *
 * @ingroup plugins
 * @author Thomas Egger
 * @date   2014-01-10
 */

class GSNIO : public IOInterface {
	public:
		GSNIO(const std::string& configfile);
		GSNIO(const GSNIO&);
		GSNIO(const Config&);

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
		void convertUnits(MeteoData& meteo) const;
		void buildStation(const std::string& id, const std::string& name, const double& lat, const double& lon,
		                  const double& alt, const double& slope_angle, const double& slope_azi, StationData &sd) const;
		bool parseMetadata(std::stringstream& ss, StationData &sd, std::string &fields, std::string &units) const;
		void readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex);
		void map_parameters(const std::string& fields, const std::string& units, MeteoData& md, std::vector<size_t>& index);
		void parse_streamElement(const std::string& line, const std::vector<size_t>& index,
		                         std::vector<MeteoData>& vecMeteo, MeteoData& tmpmeteo) const;

		void initGSNConnection();
		static size_t data_write(void* buf, size_t size, size_t nmemb, void* userp);
		bool curl_read(const std::string& url, std::ostream& os);

		const Config cfg;
		std::vector<std::string> vecStationName;
		std::map<size_t, double> multiplier, offset;
		std::string coordin, coordinparam, coordout, coordoutparam; //projection parameters
		std::string endpoint, userid, passwd; ///< Variables for endpoint configuration
		double default_timezone;
		int http_timeout; //time out for http connections
		bool gsn_debug;

		static const int http_timeout_dflt;
		static const std::string sensors_endpoint, sensors_format, null_string;
};

} //end namespace mio

#endif
