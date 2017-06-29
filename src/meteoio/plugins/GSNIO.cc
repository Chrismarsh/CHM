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
#include <meteoio/plugins/GSNIO.h>

#include <meteoio/dataClasses/Coords.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/meteoLaws/Meteoconst.h>

#include <algorithm>
#include <sstream>
#include <iostream>

#include <curl/curl.h>

using namespace std;

namespace mio {
/**
 * @page gsn GSN
 * @section gsn_format Format
 * This plugin reads meteorological data from GSN (Global Sensor Network, see <a href="http://sourceforge.net/apps/trac/gsn/"> GSN home page</a>)
 * via the RESTful web service. To compile the plugin you need to have the <a href="http://curl.haxx.se/">CURL library</a> with its headers present.
 * @subsection gsn_fields Field mapping
 * Since a lot of virtual sensors in GSN rely on parameters named with fully free text, it is necessary to map some possible names to one of MeteoIO's 
 * standard parameters. Although this plugin tries to cover as many cases as possible, it is still possible that some names would fail to be 
 * mapped to a standard parameter. In such as case, the full GSN name will be used. This means that although the data will be properly read, 
 * it will probably not be usable by a numerical model. Several solutions could then be applied:
 *     + rename the parameter in GSN in line with the names that are already recognized by this plugin;
 *     + run a small test (for example with the example "meteo_reading") in order to identify the names that are not recognized, then
 * configure MeteoIO to copy these parameter into a standard name (see \ref data_manipulations "raw data editing").
 *
 * @section gsn_units Units
 * The units of measurements are sometimes listed in the response headers, they are then parsed by the plugin and if known,
 * like <b>&deg;C</b> or <b>\%</b>, offsets and multipliers are set to convert the data to MKSA
 *
 * Otherwise the units are assumed to be the following:
 * - temperatures in celsius
 * - relative humidity in %
 * - wind speed in m/s
 * - precipitations in mm/h
 * - radiation in W/m²
 * - time is provided as a Unix timestamp, which is always in UTC
 *
 * @section gsn_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: input coordinate system (see Coords) specified in the [Input] and/or [Output] sections
 * - COORDPARAM: extra input coordinates parameters (see Coords) specified in the [Input] and/or [Output] sections
 * - GSN_URL: The URL of the RESTful web service e.g. http://montblanc.slf.ch:22001/rest
 * - GSN_USER: The username to access the service (optional)
 * - GSN_PASS: The password to authenticate the USER (optional)
 * - STATION#: station code for the given station, e. g. la_fouly_1034 (case sensitive!)
 * - GSN_TIMEOUT: timeout (in seconds) for the connection to the server (default: 60s)
 * - GSN_DEBUG: print the full requests/answers from the server when something does not work as expected
 *
 * If no STATION keys are given, the full list of ALL stations available to the user in GSN will be used!
 * This may result in a very, very long download.
 *
 * @code
 * METEO	= GSN
 * GSN_URL	= http://montblanc.slf.ch:22001/rest
 * GSN_USER	= mylogin
 * GSN_PASS	= mypasswd
 * STATION1	= wind_tunnel_meteo
 * @endcode
 *
 */

const int GSNIO::http_timeout_dflt = 60; // seconds until connect time out for libcurl
const std::string GSNIO::sensors_endpoint = "sensors";
const std::string GSNIO::sensors_format = "format=csv";
const std::string GSNIO::null_string = "null";

GSNIO::GSNIO(const std::string& configfile)
      : cfg(configfile), vecStationName(), multiplier(), offset(), coordin(),
        coordinparam(), coordout(), coordoutparam(), endpoint(), userid(), passwd(), default_timezone(1.),
        http_timeout(http_timeout_dflt), gsn_debug(false)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();

	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

GSNIO::GSNIO(const Config& cfgreader)
      : cfg(cfgreader), vecStationName(), multiplier(), offset(), coordin(),
        coordinparam(), coordout(), coordoutparam(), endpoint(), userid(), passwd(), default_timezone(1.),
        http_timeout(http_timeout_dflt), gsn_debug(false)
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	initGSNConnection();

	cfg.getValues("STATION", "INPUT", vecStationName); //reads station names into vector<string> vecStationName
}

void GSNIO::initGSNConnection() {
	curl_global_init(CURL_GLOBAL_ALL);

	cfg.getValue("GSN_TIMEOUT", "Input", http_timeout, IOUtils::nothrow);
	default_timezone = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", default_timezone, IOUtils::nothrow);

	cfg.getValue("GSN_URL", "Input", endpoint);
	if (*endpoint.rbegin() != '/') endpoint += "/";
	cerr << "[i] Using GSN URL: " << endpoint << endl;

	cfg.getValue("GSN_USER", "Input", userid, IOUtils::nothrow);
	cfg.getValue("GSN_PASS", "Input", passwd, IOUtils::nothrow);
	cfg.getValue("GSN_DEBUG", "INPUT", gsn_debug, IOUtils::nothrow);
}

void GSNIO::readStationData(const Date& date, std::vector<StationData>& vecStation)
{
	vecStation.clear();

	const size_t nrStations = vecStationName.size();
	vector<MeteoData> vecMeteo;
	for (size_t ii=0; ii<nrStations; ii++) { //loop through stations
		vecMeteo.clear();
		readData(date, date, vecMeteo, ii);
		vecStation.push_back( vecMeteo[0].meta );
	}
}

bool GSNIO::buildStation(const std::string& vs_name, const std::string& full_name, const double& lat, const double& lon,
                         const double& alt, const double& slope_angle, const double& slope_azi, StationData &sd) const
{
	Coords current_coord(coordin, coordinparam);
	if (lat==IOUtils::nodata || lon==IOUtils::nodata) 
		return false;
	
	current_coord.setLatLon(lat, lon, alt);
	const std::string name = (!full_name.empty())? full_name : vs_name;
	sd.setStationData(current_coord, vs_name, full_name);

	if (slope_angle != IOUtils::nodata) {
		if ((slope_angle == 0.) && (slope_azi == IOUtils::nodata)) {
			sd.setSlope(slope_angle, 0.); //expostion: north assumed
		} else {
			sd.setSlope(slope_angle, slope_azi);
		}
	}
	
	return true;
}

//this method is called on each station in order to parse the header and set the metadata
bool GSNIO::parseMetadata(std::stringstream& ss, StationData &sd, std::string &fields, std::string &units) const
{
	const string vsname_str("# vs_name:");
	const string altitude_str("# altitude:");
	const string longitude_str("# longitude:");
	const string latitude_str("# latitude:");
	const string slope_str("# slope:");
	const string exposition_str("# exposition:");
	const string name_str("# name:");
	const string fields_str("# fields:");
	const string units_str("# units:");

	fields.clear();
	units.clear();
	string full_name, vs_name, azi;
	double lat=IOUtils::nodata, lon=IOUtils::nodata, alt=IOUtils::nodata, slope_angle=IOUtils::nodata, slope_azi=IOUtils::nodata;
	
	string line;
	std::streamoff streampos = ss.tellg();
	while (getline(ss, line)) {
		if (line.empty() || ((line[0] != '#') && !isdigit(line[0])) ) 
			continue;

		if (isdigit(line[0])) { //reached end of metadata
			const bool status = buildStation(vs_name, full_name, lat, lon, alt, slope_angle, slope_azi, sd);
			ss.seekg(streampos, std::ios_base::beg); //point to the start of new station
			return status; //no more metadata
		}

		if (!line.compare(0, vsname_str.size(), vsname_str)) { //sensor name
			vs_name = line.substr(vsname_str.size());
			IOUtils::trim(vs_name);
		} else if (!line.compare(0, altitude_str.size(), altitude_str)) { //altitude
			IOUtils::convertString(alt, line.substr(altitude_str.size()));
		} else if (!line.compare(0, latitude_str.size(), latitude_str)) { //latitude
			IOUtils::convertString(lat, line.substr(latitude_str.size()));
		} else if (!line.compare(0, longitude_str.size(), longitude_str)) { //longitude
			IOUtils::convertString(lon, line.substr(longitude_str.size()));
		} else if (!line.compare(0, name_str.size(), name_str)) { // optional: full name
			full_name = line.substr(name_str.size());
			IOUtils::trim(full_name);
		} else if (!line.compare(0, slope_str.size(), slope_str)) { //optional: slope
			IOUtils::convertString(slope_angle, line.substr(slope_str.size()));
		} else if (!line.compare(0, exposition_str.size(), exposition_str)) { //optional: exposition
			azi = line.substr(exposition_str.size());
			if (IOUtils::isNumeric(azi)) {
				IOUtils::convertString(slope_azi, azi);
			} else {
				slope_azi = IOUtils::bearing(azi);
			}
		} else if (!line.compare(0, fields_str.size(), fields_str)) { //field
			fields = line.substr(fields_str.size());
		} else if (!line.compare(0, units_str.size(), units_str)) { //units
			units = line.substr(units_str.size()) + " "; //important when no units have been declared
		}
	}

	//This should not have happened, try to save whatever we can. Error returned by GSN will be handled afterward
	fields.clear(); //to make sure to trigger the parsing of GSN error messages
	const bool status2 = buildStation(vs_name, full_name, lat, lon, alt, slope_angle, slope_azi, sd);
	ss.seekg(streampos, std::ios_base::beg); //point to the start of new station
	return status2;
}

void GSNIO::readMeteoData(const Date& dateStart, const Date& dateEnd,
                          std::vector< std::vector<MeteoData> >& vecMeteo)
{
	const size_t nrStations = vecStationName.size();
	vecMeteo.resize(nrStations, vector<MeteoData>());
	for (size_t ii=0; ii<nrStations; ii++){ //loop through stations
		readData(dateStart, dateEnd, vecMeteo[ii], ii);
	}
}

void GSNIO::readData(const Date& dateStart, const Date& dateEnd, std::vector<MeteoData>& vecMeteo, const size_t& stationindex)
{
	const std::string station_id = vecStationName[stationindex];
	const string anon_request = sensors_endpoint + "/" + IOUtils::strToLower( station_id ) + "?" + sensors_format + "&from=" + dateStart.toString(Date::ISO)
	                            + "&to=" + dateEnd.toString(Date::ISO);
	const string auth = "&username=" + userid + "&password=" + passwd;
	const string request = (!userid.empty())? anon_request+auth : anon_request;

	stringstream ss;
	if (curl_read(request, ss)) {
		vector<size_t> index;
		string line, fields, units;

		MeteoData tmpmeteo;
		const bool meta_status = parseMetadata(ss, tmpmeteo.meta, fields, units); //read just one station

		if (units.empty() || fields.empty() || meta_status==false) {
			//when printing out a GSN error message, the # and ' ' have to be stripped from the begining -> substr(2)
			if (ss.str().find("doesn't exist in GSN!") != std::string::npos)
				throw NotFoundException(ss.str().substr(2), AT);
			if (ss.str().find("doesn't have access to the sensor") != std::string::npos)
				throw AccessException(ss.str().substr(2), AT);
			if (ss.str().find("There is no user with the provided") != std::string::npos)
				throw AccessException(ss.str().substr(2), AT);
			if (ss.str().find("The query consumed too many server resources!") != std::string::npos)
				throw IOException(ss.str().substr(2), AT);
			if (gsn_debug) {
				std::cout << "****\nRequest: " << request << "\n";
				std::cout << "Reply: " << ss.str() << "\n****\nPlease check the station name!\n";
			}
			throw InvalidFormatException("Invalid header for station " + station_id, AT);
		}
		map_parameters(fields, units, tmpmeteo, index);

		do { //parse data section, the first line should already be buffered
			if (line.empty() || (line[0] == '#') || !isdigit(line[0])) continue; //skip empty lines
			parse_streamElement(line, index, tmpmeteo);
			vecMeteo.push_back( tmpmeteo );
		} while (getline(ss, line));
		if (vecMeteo.empty()) //ie the data section was empty
			vecMeteo.push_back( tmpmeteo );
	} else {
		if (gsn_debug)
			std::cout << "****\nRequest: " << request << "\n****\n";
		throw IOException("Could not retrieve data for station " + station_id, AT);
	}
}

void GSNIO::map_parameters(const std::string& fields, const std::string& units, MeteoData& md, std::vector<size_t>& index)
{
	vector<string> field, unit;
	size_t timestamp_field = IOUtils::npos;
	multiplier.clear();
	offset.clear();

	IOUtils::readLineToVec(fields, field, ',');
	IOUtils::readLineToVec(units, unit, ',');

	if ((field.size() != unit.size()) || (field.size() < 2)) {
		throw InvalidFormatException("Fields and units are inconsistent for station " + md.meta.stationID, AT);
	}

	for (size_t ii=0; ii<field.size(); ii++) {
		const string field_name( IOUtils::strToUpper(field[ii]) );

		if (field_name == "RELATIVE_HUMIDITY" || field_name == "RH" || field_name == "AIR_HUMID" || field_name == "REL_HUMIDITY" || field_name == "RELATIVE_HUMIDITY_THYGAN" || field_name == "REL_HUMIDITY_THYGAN") {
			index.push_back(MeteoData::RH);
		} else if (field_name == "AIR_TEMPERATURE" || field_name == "TA" || field_name == "AIR_TEMP" || field_name == "AIR_TEMP_THYGAN") {
			index.push_back(MeteoData::TA);
		} else if (field_name == "WIND_DIRECTION" || field_name == "DW" || field_name == "WIND_DIRECTION_MEAN") {
			index.push_back(MeteoData::DW);
		} else if (field_name == "WIND_SPEED_MAX" || field_name == "VW_MAX") {
			index.push_back(MeteoData::VW_MAX);
		} else if (field_name == "WIND_SPEED_SCALAR_AV" || field_name == "VW" || field_name == "WIND_SPEED" || field_name == "WIND_SPEED_MEAN" || field_name == "WIND_SPEED_AV") {
			index.push_back(MeteoData::VW);
		} else if (field_name == "INCOMING_SHORTWAVE_RADIATION" || field_name == "INCOMING_SW_RADIATION" || field_name == "ISWR" || field_name == "SOLAR_RAD" || field_name == "SW_RADIATION_INCOMING") {
			index.push_back(MeteoData::ISWR);
		} else if (field_name == "INCOMING_LONGWAVE_RADIATION" || field_name == "INCOMING_LW_RADIATION" || field_name == "ILWR" || field_name == "LW_RADIATION_INCOMING") {
			index.push_back(MeteoData::ILWR);
		} else if (field_name == "OUTGOING_SHORTWAVE_RADIATION" || field_name == "RSWR") {
			index.push_back(MeteoData::RSWR);
		} else if (field_name == "OUTGOING_LONGWAVE_RADIATION" || field_name == "RLWR" || field_name == "LW_RADIATION_OUTGOING") {
			index.push_back( md.addParameter("OLWR") );
		} else if (field_name == "SNOW_HEIGHT" || field_name == "HS1") {
			index.push_back(MeteoData::HS);
		} else if (field_name == "RAIN_METER" || field_name == "PINT" || field_name == "PRECIPITATION") {
			index.push_back(MeteoData::PSUM);
		} else if (field_name == "SURFACE_TEMP" || field_name == "SURFACE_TEMPERATURE" || field_name == "TSS" || field_name == "SNOW_SURFACE_TEMPERATURE") {
			index.push_back(MeteoData::TSS);
		} else if (field_name == "SOIL_TEMP_0CM" || field_name == "SOIL_TEMP" || field_name == "SOIL_TEMPERATURE" || field_name == "TSG") {
			index.push_back(MeteoData::TSG);
		} else if (field_name == "ATM_PRESSURE" || field_name == "P") {
			index.push_back(MeteoData::P);
		} else if (field_name == "TIMESTAMP") {
			timestamp_field = ii;
			index.push_back(IOUtils::npos);
		} else if (field_name == "TIME") {
			index.push_back(IOUtils::npos);
		} else { //this is an extra parameter
			md.addParameter(field_name);
			const size_t parindex = md.getParameterIndex(field_name);
			index.push_back(parindex);

			//For the parameters unknown to MeteoIO we can store the units qualification °C, %, etc
			//and make it possible for the values to be converted to MKSA in the convertUnits procedure
			std::string name( unit[ii] );
			IOUtils::trim(name);

			if (name == "%") {
				multiplier[parindex] = 0.01;
			} else if (name.size() == 2 && (int)((unsigned char)name[0]) == 176 && name[1] == 'C') { //in °C, UTF8
				offset[parindex] = Cst::t_water_triple_pt;
			}
		}
	}

	if (timestamp_field != IOUtils::npos) { //store timestamp index at index[0]
		index[0] = timestamp_field;
	} else {
		throw InvalidFormatException("No timestamp field for station " + md.meta.stationID, AT);
	}
}

void GSNIO::parse_streamElement(const std::string& line, const std::vector<size_t>& index, MeteoData& tmpmeteo) const
{
	static vector<string> data;
	static double timestamp;
	static const size_t timestamp_index = index[0];

	const size_t size = IOUtils::readLineToVec(line, data, ',');
	if (size < 2) return; // Malformed for sure, retire gracefully, no exception thrown

	//The timestamp index is stored in index[0]
	IOUtils::convertString(timestamp, data[timestamp_index]);
	tmpmeteo.date.setUnixDate((time_t)(floor(timestamp/1000.0)));
	tmpmeteo.date.setTimeZone(default_timezone);

	size_t valid_idx=2; //index to keep track of valid parameters vs empty fields
	for (size_t jj=2; jj<size; jj++) {
		const string value = data[jj];
		if (value.empty()) continue; //skip empty values

		const size_t idx = index[valid_idx++];
		if (value != GSNIO::null_string) IOUtils::convertString(tmpmeteo(idx), value);
	}

	convertUnits(tmpmeteo);
}

void GSNIO::convertUnits(MeteoData& meteo) const
{
	//converts C to Kelvin, converts RH to [0,1]
	double& ta = meteo(MeteoData::TA);
	ta = IOUtils::C_TO_K(ta);

	double& tsg = meteo(MeteoData::TSG);
	tsg = IOUtils::C_TO_K(tsg);

	double& tss = meteo(MeteoData::TSS);
	tss = IOUtils::C_TO_K(tss);

	double& rh = meteo(MeteoData::RH);
	if (rh != IOUtils::nodata)
		rh /= 100.;

	double& hs = meteo(MeteoData::HS);
	if (hs != IOUtils::nodata)
		hs /= 100.;

	// For all parameters that have either an offset or an multiplier to bring to MKSA
	map<size_t, double>::const_iterator it;
	for (it = multiplier.begin(); it != multiplier.end(); it++) {
		double& tmp = meteo(it->first);
		if (tmp != IOUtils::nodata) tmp *= it->second;
	}

	for (it = offset.begin(); it != offset.end(); it++) {
		double& tmp = meteo(it->first);
		if (tmp != IOUtils::nodata) tmp += it->second;
	}
}

size_t GSNIO::data_write(void* buf, size_t size, size_t nmemb, void* userp)
{
	if (userp) {
		ostream& os = *static_cast<ostream*>(userp);
		const streamsize len = size * nmemb;

		if (os.write(static_cast<char*>(buf), len)) return len;
	}

	return 0;
}

bool GSNIO::curl_read(const std::string& url_query, std::ostream& os)
{
	CURLcode code(CURLE_FAILED_INIT);
	CURL* curl = curl_easy_init();

	const string url = endpoint + url_query;

	if (curl) {
		if (CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &data_write))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_FILE, &os))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_TIMEOUT, GSNIO::http_timeout))
		   && CURLE_OK == (code = curl_easy_setopt(curl, CURLOPT_URL, url.c_str())))
		{
			code = curl_easy_perform(curl);
		}
		curl_easy_cleanup(curl);
	}

	if (code!=CURLE_OK) {
		if (gsn_debug)
			std::cout << "****\nRequest: " << url_query << "\n****\n";
		std::cout << "[E] " << curl_easy_strerror(code) << "\t";
	}

	return (code==CURLE_OK);
}

} //namespace
