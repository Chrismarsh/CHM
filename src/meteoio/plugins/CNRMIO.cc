/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/plugins/CNRMIO.h>
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/MathOptim.h>
#include <meteoio/FileUtils.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/plugins/libncpp.h>

#include <cmath>
#include <cstdio>
#include <algorithm>
#include <errno.h>

using namespace std;

namespace mio {
/**
 * @page cnrm CNRM
 * @section cnrm_format Format
 * The <A HREF="http://www.cnrm.meteo.fr/">CNRM</A> has built a schema on the NetCDF format
 * to contain meteorological timeseries suitable for forcing snow models. The NetCDF (network Common Data Form) 
 * format has been created as a machine-independent format by the 
 * <A HREF="http://www.unidata.ucar.edu/">Unidata Program Center</A> in Boulder, Colorado. It is 
 * an interface for array-oriented data access and a library that provides an implementation of the interface.
 * In order to graphicaly explore the content and structure of NetCDF files, you can use the
 * <A HREF="http://www.epic.noaa.gov/java/ncBrowse/">ncBrowse</A> java software.
 *
 * @section cnrm_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] sections
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] sections
 * - METEOPATH: where to find the meteofiles as refered to here below; [Input] and [Output] sections
 * - METEOFILE: the NetCDF file which shall be used for the meteo parameter input/output (within METEOPATH); [Input] and [Output] sections
 * - UREF: height of wind measurements (in m, default 10m); [Output] section
 * - ZREF: height of air temperature measurements (in m, default 2m); [Output] section
 * - STRICTFORMAT: Whether the NetCDF file should be strictly compliant with the CNRM standard; Parameters not present
 *                 in the specification will be omitted; [Input] and [Output] section
 *
 * @section cnrm_example Example use
 * @code
 * [Input]
 * METEO     = CNRM
 * METEOPATH = ./input/meteo
 * METEOFILE = forcing.nc
 * @endcode
 * It is also recommended to use a dataGenerator on PSUM (for example, a constant generator set at 0) since
 * <A HREF="http://www.cnrm-game-meteo.fr/spip.php?article555&lang=en">Crocus</A> does not accept nodata values (and re-accumulating
 * the precipitation can still lead to nodata values).
 *
 * When using data from <A HREF="http://www.cen.ulaval.ca/nordicanad/en_index.aspx">Nordicana D</A>, it is advised to use an RhGenerator in
 * order to convert the dew point temperatures into relative humidities (and make sure the units are correct, following the steps
 * given in \ref netcdf_tricks "NetCDF tricks").
 *
 * @section cnrm_compilation Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 */

const double CNRMIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value

const std::string CNRMIO::cf_time = "time";
const std::string CNRMIO::cf_units = "units";
const std::string CNRMIO::cf_days = "days since ";
const std::string CNRMIO::cf_hours = "hours since ";
const std::string CNRMIO::cf_seconds = "seconds since ";

const std::string CNRMIO::cnrm_points = "Number_of_points";
const std::string CNRMIO::cnrm_latitude = "LAT";
const std::string CNRMIO::cnrm_longitude = "LON";
const std::string CNRMIO::cnrm_altitude = "ZS";
const std::string CNRMIO::cnrm_aspect = "aspect";
const std::string CNRMIO::cnrm_slope = "slope";
const std::string CNRMIO::cnrm_uref = "UREF";
const std::string CNRMIO::cnrm_zref = "ZREF";
const std::string CNRMIO::cnrm_ta = "Tair";
const std::string CNRMIO::cnrm_td = "d2m";
const std::string CNRMIO::cnrm_rh = "HUMREL";
const std::string CNRMIO::cnrm_vw = "Wind";
const std::string CNRMIO::cnrm_dw = "Wind_DIR";
const std::string CNRMIO::cnrm_qair = "Qair";
const std::string CNRMIO::cnrm_co2air = "CO2air";
const std::string CNRMIO::cnrm_iswr = "theorSW";
const std::string CNRMIO::cnrm_neb = "NEB";
const std::string CNRMIO::cnrm_rainf = "Rainf";
const std::string CNRMIO::cnrm_snowf = "Snowf";
const std::string CNRMIO::cnrm_swr_direct = "DIR_SWdown";
const std::string CNRMIO::cnrm_swr_diffuse = "SCA_SWdown";
const std::string CNRMIO::cnrm_p = "PSurf";
const std::string CNRMIO::cnrm_ilwr = "LWdown";
const std::string CNRMIO::cnrm_timestep = "FRC_TIME_STP";

std::map<std::string, size_t> CNRMIO::paramname;
std::map<std::string, std::string> CNRMIO::map_name;
const bool CNRMIO::__init = CNRMIO::initStaticData();

bool CNRMIO::initStaticData()
{
	//Associate unsigned int value and a string representation of a meteo parameter
	paramname[cnrm_ta] = MeteoData::TA;
	paramname[cnrm_td] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_rh] = MeteoData::RH;
	paramname[cnrm_vw] = MeteoData::VW;
	paramname[cnrm_dw] = MeteoData::DW;
	paramname[cnrm_qair] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_co2air] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_iswr] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_neb] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_rainf] = IOUtils::npos;
	paramname[cnrm_snowf] = IOUtils::npos;
	paramname[cnrm_swr_direct] = MeteoData::ISWR;
	paramname[cnrm_swr_diffuse] = IOUtils::npos;
	paramname[cnrm_p] = MeteoData::P;
	paramname[cnrm_ilwr] = MeteoData::ILWR;

	map_name["TA"] = cnrm_ta;
	map_name["RH"] = cnrm_rh;
	map_name["ILWR"] = cnrm_ilwr;
	map_name["P"] = cnrm_p;
	map_name["VW"] = cnrm_vw;
	map_name["DW"] = cnrm_dw;
	map_name["ISWR"] = cnrm_swr_direct;
	map_name["PSUM"] = cnrm_rainf;
	map_name[cnrm_co2air] = cnrm_co2air;
	map_name[cnrm_qair] = cnrm_qair;
	map_name[cnrm_neb] = cnrm_neb;

	return true;
}

CNRMIO::CNRMIO(const std::string& configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(),
                                                    in_dflt_TZ(0.), out_dflt_TZ(0.), uref(10.), zref(2.), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

CNRMIO::CNRMIO(const Config& cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(),
                                              in_dflt_TZ(0.), out_dflt_TZ(0.), uref(10.), zref(2.), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

void CNRMIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);

	cfg.getValue("STRICTFORMAT", "Input", in_strict, IOUtils::nothrow);
	cfg.getValue("STRICTFORMAT", "Output", out_strict, IOUtils::nothrow);
	
	cfg.getValue("UREF", "Output", uref, IOUtils::nothrow);
	cfg.getValue("ZREF", "Output", zref, IOUtils::nothrow);
}

void CNRMIO::readStationData(const Date&, std::vector<StationData>& vecStation)
{
	if (!vecMetaData.empty()) { // We already have meta data
		vecStation = vecMetaData;
		return;
	}

	const string path = cfg.get("METEOPATH", "Input");
	const string filename = cfg.get("METEOFILE", "Input");
	const string file_and_path = path + "/" + filename;

	if (!FileUtils::fileExists(file_and_path)) throw AccessException(file_and_path, AT); //prevent invalid filenames
	int ncid;
	ncpp::open_file(file_and_path, NC_NOWRITE, ncid);
	readMetaData(ncid, vecMetaData);
	ncpp::close_file(file_and_path, ncid);

	vecStation = vecMetaData;
}

void CNRMIO::readMetaData(const int& ncid, std::vector<StationData>& vecStation)
{
	vecStation.clear();

	int dimid;
	size_t dimlen;

	ncpp::get_dimension(ncid, cnrm_points, dimid, dimlen);
	if (dimlen == 0) return; // There are no stations

	map<string, int> map_vid;
	get_meta_data_ids(ncid, map_vid);

	double *alt = new double[dimlen];
	double *lat = new double[dimlen];
	double *lon = new double[dimlen];
	double *aspect = new double[dimlen];
	double *slope = new double[dimlen];

	ncpp::read_data(ncid, cnrm_altitude, map_vid[cnrm_altitude], alt);
	ncpp::read_data(ncid, cnrm_latitude, map_vid[cnrm_latitude], lat);
	ncpp::read_data(ncid, cnrm_longitude, map_vid[cnrm_longitude], lon);
	ncpp::read_data(ncid, cnrm_aspect, map_vid[cnrm_aspect], aspect);
	ncpp::read_data(ncid, cnrm_slope, map_vid[cnrm_slope], slope);

	//Parse to StationData objects
	Coords location(coordin, coordinparam);
	std::ostringstream ss;
	for (size_t ii=0; ii<dimlen; ii++) {
		location.setLatLon(lat[ii], lon[ii], alt[ii]);

		ss << (ii+1);
		const std::string id( ss.str() );
		ss.str("");

		ss << "Station " << (ii +1);
		const std::string name( ss.str() );
		ss.str("");

		StationData tmp(location, id, name);
		const double aspect_bearing = (aspect[ii] < 0) ? 0 : aspect[ii]; // aspect allowed to be -1 in CNRM format...
		tmp.setSlope(slope[ii], aspect_bearing);
		vecStation.push_back(tmp);
	}

	delete[] alt; delete[] lat; delete[] lon; delete[] aspect; delete[] slope;
}

void CNRMIO::readMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	vecMeteo.clear();
	const string path = cfg.get("METEOPATH", "Input");
	const string filename = cfg.get("METEOFILE", "Input");
	const string file_and_path( path + "/" + filename );

	if (!FileUtils::fileExists(file_and_path)) throw AccessException(file_and_path, AT); //prevent invalid filenames
	int ncid;
	ncpp::open_file(file_and_path, NC_NOWRITE, ncid);

	if (vecMetaData.empty()) readMetaData(ncid, vecMetaData);

	if (!vecMetaData.empty()) { //at least one station exists
		size_t index_start, index_end;
		vector<Date> vec_date;
		get_indices(ncid, dateStart, dateEnd, index_start, index_end, vec_date); //get indices for dateStart and dateEnd

		MeteoData meteo_data; //the template MeteoData object
		if ((index_start != IOUtils::npos) && (index_end != IOUtils::npos)) {
			map<string, size_t> map_parameters;
			get_parameters(ncid, map_parameters, meteo_data); //get a list of parameters present an render the template
			readData(ncid, index_start, vec_date, map_parameters, meteo_data, vecMeteo);
		}
	}

	ncpp::close_file(file_and_path, ncid);
}

void CNRMIO::readData(const int& ncid, const size_t& index_start, const std::vector<Date>& vec_date,
                        const std::map<std::string, size_t>& map_parameters, const MeteoData& meteo_data, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	const size_t number_of_stations = vecMetaData.size();
	const size_t number_of_records = vec_date.size();

	// Allocate all the MeteoData objects based on the template meteo_data
	vector<MeteoData> tmp_vec(number_of_records, meteo_data);
	for (size_t jj=0; jj<number_of_records; jj++) tmp_vec[jj].date = vec_date[jj]; //set correct date for every record

	for (size_t ii=0; ii<number_of_stations; ii++) {
		for (size_t jj=0; jj<number_of_records; jj++) tmp_vec[jj].meta = vecMetaData[ii]; //adapt meta data
		vecMeteo.push_back(tmp_vec);
	}

	// Allocate enough linear space for each parameter and read the data from NetCDF
	map<string, double*> map_data;
	for (map<string, size_t>::const_iterator it = map_parameters.begin(); it != map_parameters.end(); ++it) {
		double* data = new double[number_of_stations*number_of_records];
		const string& varname( it->first );
		map_data[varname] = data;

		int varid;
		ncpp::get_variable(ncid, varname, varid);
		ncpp::read_data_2D(ncid, varname, varid, index_start, number_of_records, number_of_stations, data);
	}

	copy_data(ncid, map_parameters, map_data, number_of_stations, number_of_records, vecMeteo);

	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); ++it) {
		delete[] it->second;
	}
}

// The copying of data into vecMeteo is a process consisting of:
// 1. A check what the relation between MeteoIO parameters and CNRM parameters present is, check map_parameters
// 2. If there is no direct association between the parameters present and the meteo_data parameters we might
//    have to deal with the parameter in a more complex way: e.g., PSUM or SWR measurements
// 3. Once we know how to deal with the parameter we loop through all stations and all parameters and copy them
//    into the appropriate places. All unit conversion have been accomplished at that point.
void CNRMIO::copy_data(const int& ncid, const std::map<std::string, size_t>& map_parameters, const std::map<std::string, double*> map_data,
                         const size_t& number_of_stations, const size_t& number_of_records, std::vector< std::vector<MeteoData> >& vecMeteo)
{
	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); ++it) {
		const string& varname( it->first );

		//find correct handling for each parameter
		bool simple_copy = false, mutiply_copy = false, psum_measurement = false, sw_measurement = false;
		double multiplier = IOUtils::nodata;
		const size_t param = map_parameters.find(varname)->second; //must exist, at this point we know it does

		if (param == IOUtils::npos) {
			if ((varname == cnrm_snowf) || (varname == cnrm_rainf)) {
				int varid;
				ncpp::get_variable(ncid, cnrm_timestep, varid);
				ncpp::read_value(ncid, cnrm_timestep, varid, multiplier);
				if (multiplier <= 0) throw InvalidArgumentException("The variable '" + cnrm_timestep + "' is invalid", AT);

				psum_measurement = true;
			} else {
				throw IOException("Don't know how to deal with parameter " + varname, AT);
			}
		} else {
			if ((varname == cnrm_swr_diffuse) || (varname == cnrm_swr_direct)) {
				sw_measurement = true;
				simple_copy = false;
			} else if (varname == cnrm_rh) {
				mutiply_copy = true;
				multiplier = 0.01;
			} else {
				simple_copy = true;
			}
		}

		// Loop through all times and all stations
		for (size_t jj=0; jj<number_of_records; jj++) {
			for (size_t ii=0; ii<number_of_stations; ii++) {
				double& value = (it->second)[jj*number_of_stations + ii];
				bool nodata = false;

				if (value == plugin_nodata) {
					nodata = true;
					value = IOUtils::nodata;
				}

				if (simple_copy) {
					vecMeteo[ii][jj](param) = value;
				} else if (mutiply_copy) {
					if (nodata) {
						vecMeteo[ii][jj](param) = value;
					} else {
						vecMeteo[ii][jj](param) = value * multiplier;
					}
				} else if (psum_measurement) { //HACK we do not handle PSUM_PH when we could use rainf / snowf
					if (!nodata) {
						double& psum = vecMeteo[ii][jj](MeteoData::PSUM);
						if (psum == IOUtils::nodata) psum = 0.0;
						psum += value * multiplier;
					}
				} else if (sw_measurement) {
					if (!nodata) {
						double& iswr = vecMeteo[ii][jj](MeteoData::ISWR);
						if (iswr == IOUtils::nodata) iswr = 0.0;
						iswr += value;
					}
				}
			}
		}
	}
}

// Go through all variables present in the NetCDF dataset that have the correct dimensions. A map called
// map_parameters will associate all parameters present with MeteoData parameters or IOUtils::npos). If
// the CNRM parameter does not have a corresponding parameter in the meteo_data object we can add a new
// parameter (e.g. cnrm_theorsw) or if the situation is more complex (e.g. rainfall is measured with two
// parameters) we deal with the situation in copy_data().
void CNRMIO::get_parameters(const int& ncid, std::map<std::string, size_t>& map_parameters, MeteoData& meteo_data)
{
	vector<string> dimensions;
	dimensions.push_back(cf_time);
	dimensions.push_back(cnrm_points);

	vector<string> parameters_present;
	ncpp::get_variables(ncid, dimensions, parameters_present);

	for (vector<string>::const_iterator it = parameters_present.begin(); it != parameters_present.end(); ++it) {
		const string name( *it );

		// Check if parameter exists in paramname, which holds strict CNRM parameters
		const map<string, size_t>::const_iterator strict_it = paramname.find(name);
		if (strict_it != paramname.end()) { // parameter is a part of the CNRM specification
			size_t index = strict_it->second;

			if (name==cnrm_qair)
				index = meteo_data.addParameter( "SH" );
			if (name==cnrm_td)
				index = meteo_data.addParameter( "TD" );
			
			if ((name == cnrm_swr_diffuse) || (name == cnrm_swr_direct) || (name == cnrm_co2air) || (name == cnrm_neb)) {
			 	index = meteo_data.addParameter( name );
			}

			map_parameters[name] = index;
		} else if (!in_strict) { // parameter will be read anyway
			size_t index = IOUtils::npos;

			if (meteo_data.param_exists(name)) {
				index = meteo_data.getParameterIndex( name );
			} else {
			 	index = meteo_data.addParameter( name );
			}

			map_parameters[name] = index;
		}
	}
}

// The CNRM format stores timestamps as doubles (either seconds or days counted from a start date)
// This method takes the dateStart and dateEnd requested and looks for the corresponding indices
// of the time variable indexStart and indexEnd.
// Furthermore the timestamps are converted to mio::Date objects and stored in vecDate
void CNRMIO::get_indices(const int& ncid, const Date& dateStart, const Date& dateEnd, size_t& indexStart, size_t& indexEnd, std::vector<Date>& vecDate)
{
	indexStart = indexEnd = IOUtils::npos;

	int varid, dimid;
	size_t dimlen;
	ncpp::get_dimension(ncid, CNRMIO::cf_time, dimid, dimlen);
	ncpp::get_variable(ncid, CNRMIO::cf_time, varid);

	// Get the units attribute and calculate the offset date
	string units_str;
	CNRMIO::TimeUnit unit_type;
	Date offset;
	ncpp::get_attribute(ncid, CNRMIO::cf_time, varid, cf_units, units_str);
	calculate_offset(units_str, unit_type, offset);

	double *time = new double[dimlen];
	ncpp::read_data(ncid, CNRMIO::cf_time, varid, time);

	// Firstly, check whether search makes any sense, that is dateStart and dateEnd overlap with the times present
	bool search = true;
	if (dimlen > 0) {
		Date time_start(offset), time_end(offset);

		double start = time[0];
		double end = time[dimlen-1];

		if (unit_type == seconds) {
			start /= 86400;
			end   /= 86400;
		}
		if (unit_type == hours) {
			start /= 24;
			end   /= 24;
		}
		time_start += Date(start, 0.0);
		time_end += Date(end, 0.0);

		if (time_start > dateEnd) search = false;
		if (time_end < dateStart) search = false;
	}

	// If search is feasible then loop through the existent timestamps and find the relevant indices
	bool start_found = false;
	if (search) {
		for (size_t ii=0; ii<dimlen; ii++) {
			if (unit_type == seconds) {
				time[ii] /= 86400;
			}
			if (unit_type == hours) {
				time[ii] /= 24;
			}

			const Date tmp_date = offset + Date(time[ii], 0.0);

			if (!start_found && (dateStart <= tmp_date && tmp_date <= dateEnd)) {
				start_found = true;
				indexStart = ii;
			} else if (start_found && (tmp_date > dateEnd)) {
				indexEnd = ii-1;
				break;
			}

			if (start_found) vecDate.push_back(tmp_date);
		}

		if (start_found && (indexEnd == IOUtils::npos)) {
			indexEnd = dimlen-1;
		}
	}

	delete[] time;
}

// The CNRM timestamps have an offset that is saved in the units attribute of
// the time variable - this method retrieves that offset
void CNRMIO::calculate_offset(const std::string& units, CNRMIO::TimeUnit& time_unit, Date& offset)
{
	string tmp(units);
	const size_t found_sec = units.find(CNRMIO::cf_seconds);
	const size_t found_hour = units.find(CNRMIO::cf_hours);
	const size_t found_day = units.find(CNRMIO::cf_days);

	if (found_sec != string::npos) {
		time_unit = seconds;
		tmp = tmp.substr(found_sec + CNRMIO::cf_seconds.size());
	} else if (found_hour != string::npos) {
		time_unit = hours;
		tmp = tmp.substr(found_hour + CNRMIO::cf_hours.size());
	} else if (found_day != string::npos) {
		time_unit = days;
		tmp = tmp.substr(found_day + CNRMIO::cf_days.size());
	} else {
		throw InvalidFormatException("Variable '"+CNRMIO::cf_time+"' has no valid attribute '" + cf_units + "'" , AT);
	}

	const bool success = IOUtils::convertString(offset, tmp, in_dflt_TZ);
	if (!success) throw InvalidFormatException("Cannot parse time: " + tmp, AT);
}

void CNRMIO::writeMeteoData(const std::vector< std::vector<MeteoData> >& vecMeteo, const std::string&)
{
	const size_t number_of_stations = vecMeteo.size();
	if (number_of_stations == 0) return; //Nothing to write

	const size_t number_of_records = vecMeteo[0].size();
	const string path = cfg.get("METEOPATH", "Output");
	const string filename = cfg.get("METEOFILE", "Output");
	const string file_and_path = path + "/" + filename;

	int ncid, did_time, vid_time, did_points;

	const bool exists = FileUtils::fileExists(file_and_path);
	if (exists) {// NOTE: file is deleted if it exists
		errno = 0;
		if ( remove(file_and_path.c_str())!=0 ) {
			std::ostringstream ss;
			ss << "Error deleting file \"" << file_and_path << "\", possible reason: " << strerror(errno);
			throw AccessException(ss.str(), AT);
		}
	}

	map<string, double*> map_data; // holds a pointer for every C array to be written
	map_data[cnrm_latitude] = new double[number_of_stations];
	map_data[cnrm_longitude] = new double[number_of_stations];
	map_data[cnrm_altitude] = new double[number_of_stations];
	map_data[cnrm_aspect] = new double[number_of_stations];
	map_data[cnrm_slope] = new double[number_of_stations];
	map_data[cnrm_uref] = new double[number_of_stations];
	map_data[cnrm_zref] = new double[number_of_stations];

	//compute data set start julian date and add Qair
	Date ref_date;
	for (size_t ii=0; ii<number_of_stations; ii++) {
		if (vecMeteo[ii].empty()) continue;
		if (ref_date.isUndef() || ref_date < vecMeteo[ii].front().date)
			ref_date = vecMeteo[ii].front().date;
	}
	ref_date.setTimeZone(0.); //so dates will be in GMT
	ref_date.rnd(3600*24, Date::DOWN); //round to the day
	
	map<string, int> varid;
	map<size_t, string> map_param_name;
	int* dates;
	get_parameters(ref_date.getJulian(), vecMeteo, map_param_name, map_data, dates);

	ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
	create_time_dimension(ref_date, ncid, did_time, vid_time);
	ncpp::add_dimension(ncid, cnrm_points, number_of_stations, did_points);
	create_meta_data(ncid, did_points, map_data, varid);
	create_parameters(ncid, did_time, did_points, number_of_records, number_of_stations, map_param_name, map_data, varid);
	ncpp::end_definitions(file_and_path, ncid);

	copy_data(number_of_stations, number_of_records, vecMeteo, map_param_name, map_data);

	ncpp::write_record(ncid, CNRMIO::cf_time, vid_time, 0, number_of_records, dates);
	for (map<string, double*>::const_iterator it = map_data.begin(); it != map_data.end(); ++it) {
		const string& varname = it->first;
		ncpp::write_data(ncid, varname, varid[varname], map_data[varname]);
		delete[] it->second;
	}

	ncpp::close_file(file_and_path, ncid);

	delete[] dates;
}

// Copy the data from the MeteoData objects into C arrays, perform all necessary
// conversions (multiplications) and set plugin_nodata values where required.
// A loop over all parameters present is performed.
void CNRMIO::copy_data(const size_t& number_of_stations, const size_t& number_of_records, const std::vector< std::vector<MeteoData> >& vecMeteo,
                         const std::map<size_t, std::string>& map_param_name, std::map<std::string, double*>& map_data_2D)
{
	for (map<size_t, string>::const_iterator it = map_param_name.begin(); it != map_param_name.end(); ++it) {
		const size_t param( it->first );
		const string varname( it->second );
		double* data = map_data_2D[varname];
		
		if (varname==cnrm_qair) { //special processing for Qair
			for (size_t ii=0; ii<number_of_stations; ++ii) {
				const double altitude = vecMeteo[ii].front().meta.position.getAltitude();
				
				for (size_t jj=0; jj<number_of_records; ++jj) {
					const double TA = vecMeteo[ii][jj](MeteoData::TA);
					const double RH = vecMeteo[ii][jj](MeteoData::RH);
					if (altitude!=IOUtils::nodata && TA!=IOUtils::nodata && RH!=IOUtils::nodata)
						data[jj*number_of_stations + ii] = Atmosphere::relToSpecHumidity(altitude, TA, RH);
					else
						data[jj*number_of_stations + ii] = plugin_nodata;
				}
			}
			continue;
		}
		
		if (varname==cnrm_rainf) {
			for (size_t ii=0; ii<number_of_stations; ++ii) {
				if (number_of_records>0) {
					map_data_2D[cnrm_rainf][ii] = plugin_nodata;
					map_data_2D[cnrm_snowf][ii] = plugin_nodata;
				}

				for (size_t jj=1; jj<number_of_records; ++jj) {
					const double acc_period = (vecMeteo[ii][jj].date.getJulian() - vecMeteo[ii][jj-1].date.getJulian()) * (24.*3600.); //in seconds
					if (acc_period==0.) continue; //this should not happen, but...
					const double psum = vecMeteo[ii][jj](MeteoData::PSUM);
					const double phase = vecMeteo[ii][jj](MeteoData::PSUM_PH);
					if (psum!=IOUtils::nodata) {
						if (phase!=IOUtils::nodata) {
							map_data_2D[cnrm_rainf][jj*number_of_stations + ii] = psum * phase / acc_period;
							map_data_2D[cnrm_snowf][jj*number_of_stations + ii] = psum * (1.-phase) / acc_period;
						} else {
							map_data_2D[cnrm_rainf][jj*number_of_stations + ii] = psum / acc_period;
							map_data_2D[cnrm_snowf][jj*number_of_stations + ii] = plugin_nodata;
						}
					} else {
							map_data_2D[cnrm_rainf][jj*number_of_stations + ii] = plugin_nodata;
							map_data_2D[cnrm_snowf][jj*number_of_stations + ii] = plugin_nodata;
					}
				}
			}
			continue;
		}
		if (varname==cnrm_snowf) continue; //this is handled by cnrm_rainf

		if (varname==cnrm_swr_diffuse) { //HACK: Crocus will redo global = dir+diff
			for (size_t ii=0; ii<number_of_stations; ++ii) {
				for (size_t jj=0; jj<number_of_records; ++jj)
					data[jj*number_of_stations + ii] = 0.;
			}
			continue;
		}

		for (size_t ii=0; ii<number_of_stations; ++ii) {
			const size_t nr_of_parameters = vecMeteo[ii].front().getNrOfParameters();
			if (param>=nr_of_parameters) { //unknown parameters are filled with nodata
				for (size_t jj=0; jj<number_of_records; ++jj)
					data[jj*number_of_stations + ii] = plugin_nodata;
				
				continue;
			}
			
			for (size_t jj=0; jj<number_of_records; ++jj) {
				const double& value = vecMeteo[ii][jj](param);

				if (value == IOUtils::nodata) {
					data[jj*number_of_stations + ii] = plugin_nodata;
				} else {
					data[jj*number_of_stations + ii] = value;
				}
			}
		}
	}
}

// Create meta data variables in the NetCDF dataset
void CNRMIO::create_meta_data(const int& ncid, const int& did, std::map<std::string, double*>& map_data_1D, std::map<std::string, int>& varid)
{
	for (map<string, double*>::const_iterator it = map_data_1D.begin(); it != map_data_1D.end(); ++it) {
		const string& varname( it->first );
		int vid;

		if (varname == cnrm_timestep) {
			ncpp::add_0D_variable(ncid, cnrm_timestep, NC_DOUBLE, vid);
		} else {
			ncpp::add_1D_variable(ncid, varname, NC_DOUBLE, did, vid);
		}
		ncpp::add_attribute(ncid, vid, "_FillValue", plugin_nodata);
		add_attributes_for_variable(ncid, vid, varname);

		varid[varname] = vid;
	}
}

// Create the parameter variables in the NetCDF dataset, allocate memory for the
// respective C arrays and store the variable ids in the varid map.
void CNRMIO::create_parameters(const int& ncid, const int& did_time, const int& did_points, const size_t& number_of_records,
                                 const size_t& number_of_stations, std::map<size_t, std::string>& map_param_name,
                                 std::map<std::string, double*>& map_data_2D, std::map<std::string, int>& varid)
{
	// At this point map_param_name holds all parameters that have values different from nodata
	for (map<size_t, string>::iterator it = map_param_name.begin(); it != map_param_name.end();) {
		bool create = false;
		string& varname( it->second );

		const map<string, string>::const_iterator it_cnrm = map_name.find(varname);
		if (it_cnrm != map_name.end()) {
			varname = it_cnrm->second; // the offical CNRM name for the parameter
			create = true;
			++it;
		} else {
			if (out_strict) {
				// ignore any parameters not defined in the CNRM standard:
				// if a parameter in map_param_name has no equivalent in the map_name map
				// it is deleted from map_param_name and henceforth ignored.
				map_param_name.erase(it++);
			} else {
				create = true;
				++it;
			}
		}

		if (create) {
			int vid;

			double* data = new double[number_of_records*number_of_stations];
			map_data_2D[varname] = data;

			ncpp::add_2D_variable(ncid, varname, NC_DOUBLE, did_time, did_points, vid);
			ncpp::add_attribute(ncid, vid, "_FillValue", plugin_nodata);
			add_attributes_for_variable(ncid, vid, varname);

			varid[varname] = vid;
		}
	}
}

// Retrieve the parameters in use (parameters, that are different from nodata
// for at least one timestamp for at least one station) and store them in
// map_param_name. map_param_name associates a MeteoData parameter index with a
// string name, that is the CNRM name for the parameter to use in the NetCDF
// file. Furthermore this method copies the meta data into the appropriate C
// arrays. The timestep interval is also calculated and added to the map_data_1D
void CNRMIO::get_parameters(const double& ref_julian, const std::vector< std::vector<MeteoData> >& vecMeteo, std::map<size_t, std::string>& map_param_name,
                              std::map<std::string, double*>& map_data_1D, int*& dates)
{
	const size_t number_of_records = vecMeteo[0].size();
	dates = new int [number_of_records];
	if (number_of_records==0) return;
	
	for (size_t ii=0; ii<number_of_records; ii++) {
		dates[ii] = static_cast<int>( Optim::round( (vecMeteo[0][ii].date.getJulian() - ref_julian) * 24.*3600. ) ); //in seconds since the start of simulation
	}

	const size_t nr_of_parameters = (!vecMeteo[0].empty())? vecMeteo[0][0].getNrOfParameters() : 0 ;
	
	vector<bool> vec_param_in_use(nr_of_parameters, false);
	vector<string> vec_param_name(nr_of_parameters, "");

	//Check consistency, dates must be existent everywhere
	bool inconsistent = false;
	for (size_t ii=0; ii<vecMeteo.size(); ++ii) {
		if (vecMeteo[ii].empty()) continue;
		
		//add metadata
		const StationData sd( vecMeteo[ii].front().meta );
		const double slope = sd.getSlopeAngle();
		const double azi = sd.getAzimuth();
		map_data_1D[cnrm_latitude][ii] = toNetcdfNodata( sd.position.getLat() );
		map_data_1D[cnrm_longitude][ii] = toNetcdfNodata( sd.position.getLon() );
		map_data_1D[cnrm_altitude][ii] = toNetcdfNodata( sd.position.getAltitude() );
		map_data_1D[cnrm_slope][ii] = (slope==IOUtils::nodata)? 0 : slope;
		map_data_1D[cnrm_aspect][ii] = (azi==IOUtils::nodata)? 0 : azi;
		map_data_1D[cnrm_uref][ii] = uref;
		map_data_1D[cnrm_zref][ii] = zref;
		
		if (number_of_records != vecMeteo[ii].size()) inconsistent = true;
		for (size_t jj=0; jj<vecMeteo[ii].size(); ++jj) {
			const MeteoData& meteo_data( vecMeteo[ii][jj] );
			
			//check time steps consistency compared to what we declared in the metadata
			const int local_date =  static_cast<int>( Optim::round( (meteo_data.date.getJulian() - ref_julian) * 24.*3600. ) );
			if (dates[jj] != local_date) inconsistent = true; //dates are ints in seconds

			//Check which parameters are in use
			for (size_t kk=0; kk<nr_of_parameters; ++kk) {
				if (!vec_param_in_use[kk]){
					if (meteo_data(kk) != IOUtils::nodata){
						vec_param_in_use[kk] = true;
						vec_param_name[kk] = meteo_data.getNameForParameter(kk);
					}
				}
			}
		}
	}

	if (inconsistent) throw IOException("Inconsistent dates in vecMeteo between different stations", AT);

	for (size_t kk=0; kk<nr_of_parameters; ++kk) {
		if (vec_param_in_use[kk])
			map_param_name[kk] = vec_param_name[kk];
	}
	map_param_name[nr_of_parameters] = cnrm_snowf; //so we can handle the precipitation splitting
	map_param_name[nr_of_parameters+1] = cnrm_qair; //add Qair
	map_param_name[nr_of_parameters+2] = cnrm_co2air; //this one needs to be there for Safran
	map_param_name[nr_of_parameters+3] = cnrm_swr_diffuse; //this one needs to be there for Safran

	double interval = (number_of_records>1)? static_cast<double>(dates[1] - dates[0]) : 0; //we consider this gives us the sampling rate
	double* timestep = new double[1];
	*timestep = interval;
	map_data_1D[cnrm_timestep] = timestep;
}

void CNRMIO::create_time_dimension(const Date& ref_date, const int& ncid, int& did_time, int& vid_time)
{
	ncpp::add_dimension(ncid, CNRMIO::cf_time, NC_UNLIMITED, did_time);
	ncpp::add_1D_variable(ncid, CNRMIO::cf_time, NC_INT, did_time, vid_time);
	ncpp::add_attribute(ncid, vid_time, "standard_name", CNRMIO::cf_time);
	ncpp::add_attribute(ncid, vid_time, "long_name", CNRMIO::cf_time);
	
	std::string ref_string( ref_date.toString(Date::ISO) );
	std::replace( ref_string.begin(), ref_string.end(), 'T', ' ');
	ncpp::add_attribute(ncid, vid_time, "units", "seconds since "+ref_string);
}

void CNRMIO::add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname)
{
	if (varname == cf_time) { //HACK this should now be deprecated
		ncpp::add_attribute(ncid, varid, "standard_name", CNRMIO::cf_time);
		ncpp::add_attribute(ncid, varid, "long_name", CNRMIO::cf_time);
		ncpp::add_attribute(ncid, varid, "units", "days since 1858-11-17 00:00:00");
	} else if (varname == CNRMIO::cnrm_timestep) {
		ncpp::add_attribute(ncid, varid, "long_name", "Forcing_Time_Step");
		ncpp::add_attribute(ncid, varid, "units", "s");
	} else if (varname == CNRMIO::cnrm_latitude) {
		ncpp::add_attribute(ncid, varid, "standard_name", "latitude");
		ncpp::add_attribute(ncid, varid, "long_name", "latitude");
		ncpp::add_attribute(ncid, varid, "units", "degrees_north");
	} else if (varname == CNRMIO::cnrm_longitude) {
		ncpp::add_attribute(ncid, varid, "standard_name", "longitude");
		ncpp::add_attribute(ncid, varid, "long_name", "longitude");
		ncpp::add_attribute(ncid, varid, "units", "degrees_east");
	} else if (varname == CNRMIO::cnrm_altitude) {
		ncpp::add_attribute(ncid, varid, "standard_name", "altitude");
		ncpp::add_attribute(ncid, varid, "long_name", "altitude");
		ncpp::add_attribute(ncid, varid, "units", "m");
		ncpp::add_attribute(ncid, varid, "positive", "up");
		ncpp::add_attribute(ncid, varid, "axis", "Z");
	} else if (varname == CNRMIO::cnrm_aspect) {
		ncpp::add_attribute(ncid, varid, "long_name", "slope aspect");
		ncpp::add_attribute(ncid, varid, "units", "degrees from north");
	} else if (varname == CNRMIO::cnrm_slope) {
		ncpp::add_attribute(ncid, varid, "long_name", "slope angle");
		ncpp::add_attribute(ncid, varid, "units", "degrees from horizontal");
	} else if (varname == CNRMIO::cnrm_uref) {
		ncpp::add_attribute(ncid, varid, "long_name", "Reference_Height_for_Wind");
		ncpp::add_attribute(ncid, varid, "units", "m");
	} else if (varname == CNRMIO::cnrm_zref) {
		ncpp::add_attribute(ncid, varid, "long_name", "Reference_Height");
		ncpp::add_attribute(ncid, varid, "units", "m");
	} else if (varname == CNRMIO::cnrm_ta) {
		ncpp::add_attribute(ncid, varid, "long_name", "Near Surface Air Temperature");
		ncpp::add_attribute(ncid, varid, "units", "K");
	} else if (varname == CNRMIO::cnrm_vw) {
		ncpp::add_attribute(ncid, varid, "long_name", "Wind Speed");
		ncpp::add_attribute(ncid, varid, "units", "m/s");
	} else if (varname == CNRMIO::cnrm_dw) {
		ncpp::add_attribute(ncid, varid, "long_name", "Wind Direction");
		ncpp::add_attribute(ncid, varid, "units", "deg");
	} else if (varname == CNRMIO::cnrm_iswr) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Incident total Shortwave radiation");
		ncpp::add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == CNRMIO::cnrm_swr_direct) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Incident Direct Shortwave Radiation");
		ncpp::add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == CNRMIO::cnrm_swr_diffuse) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Incident Diffuse Shortwave Radiation");
		ncpp::add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == CNRMIO::cnrm_rainf) {
		ncpp::add_attribute(ncid, varid, "long_name", "Rainfall Rate");
		ncpp::add_attribute(ncid, varid, "units", "kg/m2/s");
	} else if (varname == CNRMIO::cnrm_snowf) {
		ncpp::add_attribute(ncid, varid, "long_name", "Snowfall Rate");
		ncpp::add_attribute(ncid, varid, "units", "kg/m2/s");
	} else if (varname == CNRMIO::cnrm_rh) {
		ncpp::add_attribute(ncid, varid, "long_name", "Relative Humidity");
		ncpp::add_attribute(ncid, varid, "units", "%");
	} else if (varname == CNRMIO::cnrm_qair) {
		ncpp::add_attribute(ncid, varid, "long_name", "Near Surface Specific Humidity");
		ncpp::add_attribute(ncid, varid, "units", "Kg/Kg");
	} else if (varname == CNRMIO::cnrm_ilwr) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Incident Longwave Radiation");
		ncpp::add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == CNRMIO::cnrm_p) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Pressure");
		ncpp::add_attribute(ncid, varid, "units", "Pa");
	} else if (varname == CNRMIO::cnrm_co2air) {
		ncpp::add_attribute(ncid, varid, "long_name", "Near Surface CO2 Concentration");
		ncpp::add_attribute(ncid, varid, "units", "kg/m3");
	}
}

void CNRMIO::get_meta_data_ids(const int& ncid, std::map<std::string, int>& map_vid)
{
	const string names[] = {cnrm_altitude, cnrm_latitude, cnrm_longitude, cnrm_aspect, cnrm_slope};
	vector<string> varname(names, names + sizeof(names) / sizeof(names[0]));

	vector<string> dimensions;
	dimensions.push_back(cnrm_points); // All variables have to have the dimension cnrm_points

	for (vector<string>::const_iterator it = varname.begin(); it != varname.end(); ++it) {
		int varid;
		const string& name( *it );

		ncpp::get_variable(ncid, name, varid);
		if (!ncpp::check_dimensions(ncid, name, varid, dimensions))
			throw IOException("Variable '" + name  + "' fails dimension check", AT);

		map_vid[name] = varid;
	}
}

double CNRMIO::toNetcdfNodata(const double& value) const
{
		if (value==IOUtils::nodata) return plugin_nodata;
		else return value;
}

} //namespace
