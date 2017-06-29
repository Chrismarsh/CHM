/***********************************************************************************/
/*  Copyright 2009 SLF                                                                                                                                */
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
#ifndef GSNIO_H
#define GSNIO_H

#include <meteoio/IOInterface.h>

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

		virtual void readStationData(const Date& date, std::vector<StationData>& vecStation);
		virtual void readMeteoData(const Date& dateStart, const Date& dateEnd,
		                           std::vector< std::vector<MeteoData> >& vecMeteo);

	private:
		void convertUnits(MeteoData& meteo) const;
		bool buildStation(const std::string& vs_name, const std::string& full_name, const double& lat, const double& lon,
		                  const double& alt, const double& slope_angle, const double& slope_azi, StationData &sd) const;
		bool parseMetadata(std::stringstream& ss, StationData &sd, std::string &fields, std::string &units) const;
		void readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex);
		void map_parameters(const std::string& fields, const std::string& units, MeteoData& md, std::vector<size_t>& index);
		void parse_streamElement(const std::string& line, const std::vector<size_t>& index, MeteoData& tmpmeteo) const;

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
