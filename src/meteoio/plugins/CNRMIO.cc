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
#include "CNRMIO.h"
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/Timer.h>
#include <meteoio/MathOptim.h>
#include <meteoio/plugins/libncpp.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

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
 * @section cnrm_units Units
 *
 *
 * @section cnrm_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] sections
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] sections
 * - DEMFILE: The filename of the file containing the DEM; [Input] section
 * - DEMVAR: The variable name of the DEM within the DEMFILE; [Input] section
 * - METEOPATH: where to find the meteofiles as refered to here below; [Input] and [Output] sections
 * - METEOFILE: the NetCDF file which shall be used for the meteo parameter input/output; [Input] and [Output] sections
 * - GRID2DFILE: the NetCDF file which shall be used for gridded input/output; [Input] and [Output] sections
 * - STRICTFORMAT: Whether the NetCDF file should be strictly compliant with the CNRM standard; Parameters not present
 *                 in the specification will be omitted; [Input] and [Output] section
 *
 * @section cnrm_example Example use
 * @code
 * [Input]
 * METEO     = CNRM
 * METEOFILE = ./input/meteo/forcing.nc
 * @endcode
 *
 * @section cnrm_compilation Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 */

const double CNRMIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value
const double CNRMIO::epsilon = 1.0e-10; //when comparing timestamps

const std::string CNRMIO::cf_time = "time";
const std::string CNRMIO::cf_units = "units";
const std::string CNRMIO::cf_days = "days since ";
const std::string CNRMIO::cf_hours = "hours since ";
const std::string CNRMIO::cf_seconds = "seconds since ";
const std::string CNRMIO::cf_latitude = "lat";
const std::string CNRMIO::cf_longitude = "lon";
const std::string CNRMIO::cf_altitude = "z";
const std::string CNRMIO::cf_ta = "temperature";
const std::string CNRMIO::cf_rh = "humidity";
const std::string CNRMIO::cf_p = "pressure";

const std::string CNRMIO::cnrm_points = "Number_of_points";
const std::string CNRMIO::cnrm_latitude = "LAT";
const std::string CNRMIO::cnrm_longitude = "LON";
const std::string CNRMIO::cnrm_altitude = "ZS";
const std::string CNRMIO::cnrm_aspect = "aspect";
const std::string CNRMIO::cnrm_slope = "slope";
const std::string CNRMIO::cnrm_ta = "Tair";
const std::string CNRMIO::cnrm_rh = "HUMREL";
const std::string CNRMIO::cnrm_vw = "Wind";
const std::string CNRMIO::cnrm_dw = "Wind_DIR";
const std::string CNRMIO::cnrm_qair = "Qair";
const std::string CNRMIO::cnrm_co2air = "CO2air";
const std::string CNRMIO::cnrm_theorsw = "theorSW";
const std::string CNRMIO::cnrm_neb = "NEB";
const std::string CNRMIO::cnrm_psum = "Rainf";
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
	paramname[cnrm_qair] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_co2air] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_neb] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_theorsw] = IOUtils::npos; // not a standard MeteoIO parameter
	paramname[cnrm_rh] = MeteoData::RH;
	paramname[cnrm_vw] = MeteoData::VW;
	paramname[cnrm_dw] = MeteoData::DW;
	paramname[cnrm_psum] = IOUtils::npos;
	paramname[cnrm_snowf] = IOUtils::npos;
	paramname[cnrm_swr_direct] = IOUtils::npos;
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
	map_name["PSUM"] = cnrm_psum;
	map_name[cnrm_co2air] = cnrm_co2air;
	map_name[cnrm_qair] = cnrm_qair;
	map_name[cnrm_theorsw] = cnrm_theorsw;
	map_name[cnrm_neb] = cnrm_neb;

	return true;
}

CNRMIO::CNRMIO(const std::string& configfile) : cfg(configfile), coordin(), coordinparam(), coordout(), coordoutparam(),
                                                    in_dflt_TZ(0.), out_dflt_TZ(0.), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

CNRMIO::CNRMIO(const Config& cfgreader) : cfg(cfgreader), coordin(), coordinparam(), coordout(), coordoutparam(),
                                              in_dflt_TZ(0.), out_dflt_TZ(0.), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

CNRMIO::~CNRMIO() throw() {}

void CNRMIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);

	cfg.getValue("STRICTFORMAT", "Input", in_strict, IOUtils::nothrow);
	cfg.getValue("STRICTFORMAT", "Output", out_strict, IOUtils::nothrow);
}

void CNRMIO::read2DGrid(Grid2DObject& grid_out, const std::string& arguments)
{
	vector<string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() == 2) {
		read2DGrid_internal(grid_out, vec_argument[0], vec_argument[1]);
	} else {
		throw InvalidArgumentException("The format for the arguments to CNRMIO::read2DGrid is filename:varname", AT);
	}
}

void CNRMIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const string filename = cfg.get("GRID2DFILE", "Input");
	const string varname = get_varname(parameter);
	read2DGrid_internal(grid_out, filename, varname, date);
}

void CNRMIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& filename, const std::string& varname, const Date& date)
{
	const bool is_record = (date != Date());
	size_t lat_index = 0, lon_index = 1;

	int ncid, varid;
	vector<int> dimid, dim_varid;
	vector<string> dimname;
	vector<size_t> dimlen;

	if (!IOUtils::fileExists(filename)) throw FileAccessException(filename, AT); //prevent invalid filenames
	ncpp::open_file(filename, NC_NOWRITE, ncid);
	ncpp::get_variable(ncid, varname, varid);
	ncpp::get_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);

	if (is_record) { // In case we're reading a record the first index is always the record index
		lat_index = 1;
		lon_index = 2;

		if (dimid.size()!=3 || dimlen[0]<1 || dimlen[lat_index]<2 || dimlen[lon_index]<2)
			throw IOException("Variable '" + varname + "' may only have three dimensions, all have to at least have length 1", AT);
	} else if (dimid.size()==3 && dimlen[0]==1) { //in case the variable is associated with a 1 element time dimension
		lat_index = 1;
		lon_index = 2;

		if (dimlen[lat_index]<2 || dimlen[lon_index]<2)
			throw IOException("All dimensions for variable '" + varname + "' have to at least have length 1", AT);
	} else if (dimid.size()!=2 || dimlen[lat_index]<2 || dimlen[lon_index]<2) {
		throw IOException("Variable '" + varname + "' may only have two dimensions and both have to have length >1", AT);
	}

	double *lat = new double[dimlen[lat_index]];
	double *lon = new double[dimlen[lon_index]];
	double *grid = new double[dimlen[lat_index]*dimlen[lon_index]];

	ncpp::read_data(ncid, dimname[lat_index], dim_varid[lat_index], lat);
	ncpp::read_data(ncid, dimname[lon_index], dim_varid[lon_index], lon);

	if (is_record) {
		const size_t pos = ncpp::find_record(ncid, CNRMIO::cf_time, dimid[0], date.getModifiedJulianDate());
		if (pos == IOUtils::npos)
			throw IOException("No record for date " + date.toString(Date::ISO), AT);

		ncpp::read_data(ncid, varname, varid, pos, dimlen[lat_index], dimlen[lon_index], grid);
	} else {
		ncpp::read_data(ncid, varname, varid, grid);
	}

	double missing_value=plugin_nodata;
	if (ncpp::check_attribute(ncid, varid, "missing_value"))
		ncpp::get_attribute(ncid, varname, varid, "missing_value", missing_value);

	ncpp::copy_grid(coordin, coordinparam, dimlen[lat_index], dimlen[lon_index], lat, lon, grid, missing_value, grid_out);

	//handle data packing if necessary
	if (ncpp::check_attribute(ncid, varid, "scale_factor")) {
		double scale_factor=1.;
		ncpp::get_attribute(ncid, varname, varid, "scale_factor", scale_factor);
		grid_out.grid2D *= scale_factor;
	}
	if (ncpp::check_attribute(ncid, varid, "add_offset")) {
		double add_offset=0.;
		ncpp::get_attribute(ncid, varname, varid, "add_offset", add_offset);
		grid_out.grid2D += add_offset;
	}

	ncpp::close_file(filename, ncid);
	delete[] lat; delete[] lon; delete[] grid;
}

void CNRMIO::readDEM(DEMObject& dem_out)
{
	const string filename = cfg.get("DEMFILE", "Input");
	const string varname = cfg.get("DEMVAR", "Input");
	read2DGrid_internal(dem_out, filename, varname);
}

void CNRMIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CNRMIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
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

	int ncid;
	if (!IOUtils::fileExists(file_and_path)) throw FileAccessException(file_and_path, AT); //prevent invalid filenames
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
	map<string, int> map_vid;

	ncpp::get_dimension(ncid, cnrm_points, dimid, dimlen);
	if (dimlen == 0) return; // There are no stations

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
	ostringstream ss;
	for (size_t ii=0; ii<dimlen; ii++) {
		location.setLatLon(lat[ii], lon[ii], alt[ii]);

		ss << (ii+1);
		const string id( ss.str() );
		ss.str("");

		ss << "Station " << (ii +1);
		const string name( ss.str() );
		ss.str("");

		StationData tmp(location, id, name);
		const double aspect_bearing = (aspect[ii] < 0) ? 0 : aspect[ii]; // aspect allowed to be -1 in CNRM format...
		tmp.setSlope(slope[ii], aspect_bearing);
		vecStation.push_back(tmp);
	}

	delete[] alt; delete[] lat; delete[] lon; delete[] aspect; delete[] slope;
}

void CNRMIO::readMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< std::vector<MeteoData> >& vecMeteo, const size_t&)
{
	vecMeteo.clear();
	const string path = cfg.get("METEOPATH", "Input");
	const string filename = cfg.get("METEOFILE", "Input");
	const string file_and_path = path + "/" + filename;

	int ncid;
	if (!IOUtils::fileExists(file_and_path)) throw FileAccessException(file_and_path, AT); //prevent invalid filenames
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
		const string& varname = it->first;

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
		const string& varname = it->first;

		//find correct handling for each parameter
		bool simple_copy = false, mutiply_copy = false, psum_measurement = false, sw_measurement = false;
		double multiplier = IOUtils::nodata;
		const size_t param = map_parameters.find(varname)->second; //must exist, at this point we know it does

		if (param == IOUtils::npos) {
			if ((varname == cnrm_snowf) || (varname == cnrm_psum)) {
				int varid;
				ncpp::get_variable(ncid, cnrm_timestep, varid);
				ncpp::read_value(ncid, cnrm_timestep, varid, multiplier);

				if (multiplier <= 0) throw InvalidArgumentException("The variable '" + cnrm_timestep + "' is invalid", AT);

				psum_measurement = true;
			} else if ((varname == cnrm_swr_diffuse) || (varname == cnrm_swr_direct)) {
				sw_measurement = true;
			} else {
				throw IOException("Don't know how to deal with parameter " + varname, AT);
			}
		} else {
			if (varname == cnrm_rh) {
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
				} else if (psum_measurement) {
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
		const string& name = *it;
		//cout << "Found parameter: " << name << endl;

		// Check if parameter exists in paramname, which holds strict CNRM parameters
		const map<string, size_t>::const_iterator strict_it = paramname.find(name);
		if (strict_it != paramname.end()) { // parameter is a part of the CNRM specification
			size_t index = strict_it->second;

			if ((name == cnrm_theorsw) || (name == cnrm_qair) || (name == cnrm_co2air) || (name == cnrm_neb)) {
			 	index = meteo_data.addParameter(name);
			}

			map_parameters[name] = index;
		} else if (!in_strict) { // parameter will be read anyway
			size_t index = IOUtils::npos;

			if (meteo_data.param_exists(name)) {
				index = meteo_data.getParameterIndex(name);
			} else {
			 	index = meteo_data.addParameter(name);
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
	bool create_time = false, create_points = false, create_locations = false, create_variables = false;

	const bool exists = IOUtils::fileExists(file_and_path);
	if (exists) remove(file_and_path.c_str()); // NOTE: file is deleted if it exists

	double* dates;
	map<string, double*> map_data; // holds a pointer for every C array to be written
	map_data[cnrm_latitude] = new double[number_of_stations];
	map_data[cnrm_longitude] = new double[number_of_stations];
	map_data[cnrm_altitude] = new double[number_of_stations];
	map_data[cnrm_aspect] = new double[number_of_stations];
	map_data[cnrm_slope] = new double[number_of_stations];

	map<string, int> varid;
	map<size_t, string> map_param_name;

	get_parameters(vecMeteo, map_param_name, map_data, dates);

	ncpp::create_file(file_and_path, NC_CLASSIC_MODEL, ncid);
	create_time = create_points = create_locations = create_variables = true;

	if (create_time) create_time_dimension(ncid, did_time, vid_time);
	if (create_points) ncpp::add_dimension(ncid, cnrm_points, number_of_stations, did_points);
	if (create_locations) create_meta_data(ncid, did_points, map_data, varid);
	if (create_variables) create_parameters(ncid, did_time, did_points, number_of_records, number_of_stations, map_param_name, map_data, varid);

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
		const size_t param = it->first;
		const string varname = it->second;

		bool simple_copy = false, multiply_copy = false;
		double multiplier = IOUtils::nodata;

		double* data = map_data_2D[varname];

		if (param == MeteoData::RH) {
			multiplier = 100.;
			multiply_copy = true;
		} else if (param == MeteoData::PSUM) {
			multiply_copy = true;
			multiplier = 1./3600.;
		} else {
			simple_copy = true;
		}

		for (size_t ii=0; ii<number_of_stations; ++ii) {
			for (size_t jj=0; jj<number_of_records; ++jj) {
				const double& value = vecMeteo[ii][jj](param);

				if (value == IOUtils::nodata) {
					data[jj*number_of_stations + ii] = plugin_nodata;
				} else if (simple_copy) {
					data[jj*number_of_stations + ii] = value;
				} else if (multiply_copy) {
					data[jj*number_of_stations + ii] = value * multiplier;
				}
			}
		}
	}
}

// Create meta data variables in the NetCDF dataset
void CNRMIO::create_meta_data(const int& ncid, const int& did, std::map<std::string, double*>& map_data_1D, std::map<std::string, int>& varid)
{
	for (map<string, double*>::const_iterator it = map_data_1D.begin(); it != map_data_1D.end(); ++it) {
		int vid;
		const string& varname = it->first;

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
		string& varname = it->second;

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
void CNRMIO::get_parameters(const std::vector< std::vector<MeteoData> >& vecMeteo, std::map<size_t, std::string>& map_param_name,
                              std::map<std::string, double*>& map_data_1D, double*& dates)
{
	const size_t number_of_records = vecMeteo[0].size();
	dates = new double[number_of_records];

	double interval = 0;
	for (size_t ii=0; ii<number_of_records; ii++) {
		dates[ii] = vecMeteo[0][ii].date.getModifiedJulianDate();
		if (ii == 1) interval = static_cast<double>( Optim::round((dates[ii] - dates[ii-1]) * 86400.) );
	}

	const size_t nr_of_parameters = (!vecMeteo[0].empty())? vecMeteo[0][0].getNrOfParameters() : 0 ;
	vector<bool> vec_param_in_use(nr_of_parameters, false);
	vector<string> vec_param_name(nr_of_parameters, "");

	//Check consistency, dates must be existent everywhere
	bool inconsistent = false;
	for (size_t ii=0; ii<vecMeteo.size(); ++ii) {
		if (number_of_records != vecMeteo[ii].size()) inconsistent = true;
		for (size_t jj=0; jj<vecMeteo[ii].size(); ++jj) {
			const MeteoData& meteo_data = vecMeteo[ii][jj];

			if (!IOUtils::checkEpsilonEquality(dates[jj], meteo_data.date.getModifiedJulianDate(), CNRMIO::epsilon))
				inconsistent = true;

			if (jj == 0) {
				map_data_1D[cnrm_latitude][ii] = meteo_data.meta.position.getLat();
				map_data_1D[cnrm_longitude][ii] = meteo_data.meta.position.getLon();
				map_data_1D[cnrm_altitude][ii] = meteo_data.meta.position.getAltitude();
				map_data_1D[cnrm_slope][ii] = meteo_data.meta.getSlopeAngle();
				map_data_1D[cnrm_aspect][ii] = meteo_data.meta.getAzimuth();
			}

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


	double* timestep = new double[1];
	*timestep = interval;
	map_data_1D[cnrm_timestep] = timestep;
}

void CNRMIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void CNRMIO::write2DGrid(const Grid2DObject& grid_in, const std::string& arguments)
{
	// arguments is a string of the format filname:varname
	vector<string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() != 2)
		throw InvalidArgumentException("The format for the arguments to CNRMIO::write2DGrid is filename:varname", AT);

	write2DGrid_internal(grid_in, vec_argument[0], vec_argument[1]);
}

void CNRMIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const string filename = cfg.get("GRID2DFILE", "Output");
	const string varname = get_varname(parameter);

	write2DGrid_internal(grid_in, filename, varname, date);
}

void CNRMIO::write2DGrid_internal(const Grid2DObject& grid_in, const std::string& filename, const std::string& varname, const Date& date)
{
	const bool is_record = (date != Date());
	const bool exists = IOUtils::fileExists(filename);

	double *lat_array = new double[grid_in.getNy()];
	double *lon_array = new double[grid_in.getNx()];
	double *data = new double[grid_in.getNy() * grid_in.getNx()];

	ncpp::calculate_dimensions(grid_in, lat_array, lon_array);
	ncpp::fill_grid_data(grid_in, data);

	int ncid, did_lat, did_lon, did_time, vid_lat, vid_lon, vid_var, vid_time;
	bool create_dimensions(false), create_variable(false), create_time(false);

	if (exists) {
		ncpp::open_file(filename, NC_WRITE, ncid);

		//check of lat/lon are defined and consistent
		if (ncpp::check_dim_var(ncid, cf_latitude) && ncpp::check_dim_var(ncid, cf_longitude)) {
			check_consistency(ncid, grid_in, lat_array, lon_array, did_lat, did_lon, vid_lat, vid_lon);
		} else {
			create_dimensions = true;
		}

		if (is_record) {
			//check if a time dimension/variable already exists
			if (ncpp::check_dim_var(ncid, CNRMIO::cf_time)) {
				ncpp::get_dimension(ncid, CNRMIO::cf_time, did_time);
				ncpp::get_variable(ncid, CNRMIO::cf_time, vid_time);
			} else {
				create_time = true;
			}
		}

		if (ncpp::check_variable(ncid, varname)) { // variable exists
			ncpp::get_variable(ncid, varname, vid_var);

			vector<int> dimid, dim_varid;
			vector<string> dimname;
			vector<size_t> dimlen;

			ncpp::get_dimension(ncid, varname, vid_var, dimid, dim_varid, dimname, dimlen);

			if (is_record) {
				if ((dimname.size() != 3) || (dimname[0] != cf_time) || (dimname[1] != cf_latitude) || (dimname[2] != cf_longitude) || (dimlen[1]!=grid_in.getNy()) || (dimlen[2]!=grid_in.getNx()))
					throw IOException("Variable '" + varname  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
			} else {
				if ((dimname[0] != cf_latitude) || (dimname[1] != cf_longitude) || (dimlen[0]!=grid_in.getNy()) || (dimlen[1]!=grid_in.getNx()))
					throw IOException("Variable '" + varname  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
			}
		} else {
			create_variable = true;
		}

		ncpp::start_definitions(filename, ncid);
	} else {
		if (!IOUtils::validFileAndPath(filename)) throw InvalidFileNameException(filename,AT);
		ncpp::create_file(filename, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.3");

		create_variable = create_dimensions = true;

		if (is_record) create_time = true;
	}

	if (create_dimensions) create_latlon_dimensions(ncid, grid_in, did_lat, did_lon, vid_lat, vid_lon);
	if (create_time) create_time_dimension(ncid, did_time, vid_time);

	if (is_record && create_variable) {
		ncpp::add_3D_variable(ncid, varname, NC_DOUBLE, did_time, did_lat, did_lon, vid_var);
		add_attributes_for_variable(ncid, vid_var, varname);
	} else if (create_variable) {
		ncpp::add_2D_variable(ncid, varname, NC_DOUBLE, did_lat, did_lon, vid_var);
		add_attributes_for_variable(ncid, vid_var, varname);
	}

	ncpp::end_definitions(filename, ncid);

	if (create_dimensions) {
		ncpp::write_data(ncid, cf_latitude, vid_lat, lat_array);
		ncpp::write_data(ncid, cf_longitude, vid_lon, lon_array);
	}

	if (is_record) {
		size_t pos_start = ncpp::add_record(ncid, CNRMIO::cf_time, vid_time, date.getModifiedJulianDate());
		ncpp::write_data(ncid, varname, vid_var, grid_in.getNy(), grid_in.getNx(), pos_start, data);
	} else {
		ncpp::write_data(ncid, varname, vid_var, data);
	}

	ncpp::close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

void CNRMIO::create_latlon_dimensions(const int& ncid, const Grid2DObject& grid_in, int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	ncpp::add_dimension(ncid, cf_latitude, grid_in.getNy(), did_lat);
	ncpp::add_1D_variable(ncid, cf_latitude, NC_DOUBLE, did_lat, vid_lat);
	add_attributes_for_variable(ncid, vid_lat, cf_latitude);

	ncpp::add_dimension(ncid, cf_longitude, grid_in.getNx(), did_lon);
	ncpp::add_1D_variable(ncid, cf_longitude, NC_DOUBLE, did_lon, vid_lon);
	add_attributes_for_variable(ncid, vid_lon, cf_longitude);
}

void CNRMIO::create_time_dimension(const int& ncid, int& did_time, int& vid_time)
{
	ncpp::add_dimension(ncid, CNRMIO::cf_time, NC_UNLIMITED, did_time);
	ncpp::add_1D_variable(ncid, CNRMIO::cf_time, NC_DOUBLE, did_time, vid_time); // julian day
	add_attributes_for_variable(ncid, vid_time, CNRMIO::cf_time);
}

// When reading or writing gridded variables we should have a consistent naming
// scheme: http://cfconventions.org/1.6.html
std::string CNRMIO::get_varname(const MeteoGrids::Parameters& parameter)
{
	string varname = MeteoGrids::getParameterName(parameter);

	if (parameter == MeteoGrids::TA) varname = cnrm_ta;
	else if (parameter == MeteoGrids::RH) varname = cnrm_rh;
	else if (parameter == MeteoGrids::DEM) varname = cnrm_altitude;
	else if (parameter == MeteoGrids::P) varname = cnrm_p;
	else if (parameter == MeteoGrids::VW) varname = cnrm_vw;
	else if (parameter == MeteoGrids::DW) varname = cnrm_dw;
	else if (parameter == MeteoGrids::ILWR) varname = cnrm_ilwr;
	else if (parameter == MeteoGrids::PSUM) varname = cnrm_psum; //HACK this should add snowf!
	else if (parameter == MeteoGrids::SLOPE) varname = cnrm_slope;
	else if (parameter == MeteoGrids::AZI) varname = cnrm_aspect;
	//HACK: iswr=dir+diff
	//HACK: U, V from vw, dw

	return varname;
}

void CNRMIO::add_attributes_for_variable(const int& ncid, const int& varid, const std::string& varname)
{
	if (varname == cf_latitude) {
		ncpp::add_attribute(ncid, varid, "standard_name", "latitude");
		ncpp::add_attribute(ncid, varid, "long_name", "latitude");
		ncpp::add_attribute(ncid, varid, "units", "degrees_north");
	} else if (varname == cf_longitude) {
		ncpp::add_attribute(ncid, varid, "standard_name", "longitude");
		ncpp::add_attribute(ncid, varid, "long_name", "longitude");
		ncpp::add_attribute(ncid, varid, "units", "degrees_east");
	} else if (varname == cf_altitude) {
		ncpp::add_attribute(ncid, varid, "standard_name", "altitude");
		ncpp::add_attribute(ncid, varid, "long_name", "height above mean sea level");
		ncpp::add_attribute(ncid, varid, "units", "m");
		ncpp::add_attribute(ncid, varid, "positive", "up");
		ncpp::add_attribute(ncid, varid, "axis", "Z");
	} else if (varname == cf_p) {
		ncpp::add_attribute(ncid, varid, "standard_name", "air_pressure");
		ncpp::add_attribute(ncid, varid, "long_name", "near surface air pressure");
		ncpp::add_attribute(ncid, varid, "units", "Pa");
	} else if (varname == cf_ta) {
		ncpp::add_attribute(ncid, varid, "standard_name", "air_temperature");
		ncpp::add_attribute(ncid, varid, "long_name", "near surface air temperature");
		ncpp::add_attribute(ncid, varid, "units", "K");
	} else if (varname == cf_rh) {
		ncpp::add_attribute(ncid, varid, "standard_name", "relative humidity");
		ncpp::add_attribute(ncid, varid, "long_name", "relative humidity");
		ncpp::add_attribute(ncid, varid, "units", "fraction");
	} else if (varname == cf_time) {
		ncpp::add_attribute(ncid, varid, "standard_name", CNRMIO::cf_time);
		ncpp::add_attribute(ncid, varid, "long_name", CNRMIO::cf_time);
		ncpp::add_attribute(ncid, varid, "units", "days since 1858-11-17 00:00:00");
	} else if (varname == CNRMIO::cnrm_altitude) {
		ncpp::add_attribute(ncid, varid, "long_name", "altitude");
		ncpp::add_attribute(ncid, varid, "units", "m");
	} else if (varname == CNRMIO::cnrm_aspect) {
		ncpp::add_attribute(ncid, varid, "long_name", "slope aspect");
		ncpp::add_attribute(ncid, varid, "units", "degrees from north");
	} else if (varname == CNRMIO::cnrm_slope) {
		ncpp::add_attribute(ncid, varid, "long_name", "slope angle");
		ncpp::add_attribute(ncid, varid, "units", "degrees from horizontal");
	} else if (varname == CNRMIO::cnrm_latitude) {
		ncpp::add_attribute(ncid, varid, "long_name", "latitude");
		ncpp::add_attribute(ncid, varid, "units", "degrees_north");
	} else if (varname == CNRMIO::cnrm_longitude) {
		ncpp::add_attribute(ncid, varid, "long_name", "longitude");
		ncpp::add_attribute(ncid, varid, "units", "degrees_east");
	} else if (varname == CNRMIO::cnrm_ta) {
		ncpp::add_attribute(ncid, varid, "long_name", "Near Surface Air Temperature");
		ncpp::add_attribute(ncid, varid, "units", "K");
	} else if (varname == CNRMIO::cnrm_timestep) {
		ncpp::add_attribute(ncid, varid, "long_name", "Forcing_Time_Step");
		ncpp::add_attribute(ncid, varid, "units", "s");
	} else if (varname == CNRMIO::cnrm_vw) {
		ncpp::add_attribute(ncid, varid, "long_name", "Wind Speed");
		ncpp::add_attribute(ncid, varid, "units", "m/s");
	} else if (varname == CNRMIO::cnrm_dw) {
		ncpp::add_attribute(ncid, varid, "long_name", "Wind Direction");
		ncpp::add_attribute(ncid, varid, "units", "deg");
	} else if (varname == CNRMIO::cnrm_swr_direct) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Incident Direct Shortwave Radiation");
		ncpp::add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == CNRMIO::cnrm_psum) {
		ncpp::add_attribute(ncid, varid, "long_name", "Rainfall Rate");
		ncpp::add_attribute(ncid, varid, "units", "kg/m2/s");
	} else if (varname == CNRMIO::cnrm_rh) {
		ncpp::add_attribute(ncid, varid, "long_name", "Relative Humidity");
		ncpp::add_attribute(ncid, varid, "units", "%");
	} else if (varname == CNRMIO::cnrm_ilwr) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Incident Longwave Radiation");
		ncpp::add_attribute(ncid, varid, "units", "W/m2");
	} else if (varname == CNRMIO::cnrm_p) {
		ncpp::add_attribute(ncid, varid, "long_name", "Surface Pressure");
		ncpp::add_attribute(ncid, varid, "units", "Pa");
	}
}

void CNRMIO::check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
                                 int& did_lat, int& did_lon, int& vid_lat, int& vid_lon)
{
	size_t latlen, lonlen;

	ncpp::get_dimension(ncid, cf_latitude, did_lat, latlen);
	ncpp::get_dimension(ncid, cf_longitude, did_lon, lonlen);

	ncpp::get_variable(ncid, cf_latitude, vid_lat);
	ncpp::get_variable(ncid, cf_longitude, vid_lon);

	if ((latlen != grid.getNy()) || (lonlen != grid.getNx()))
		throw IOException("Error while writing grid - grid size and lat/lon coordinates are inconsistent", AT);

	double *lat = new double[grid.getNy()];
	double *lon = new double[grid.getNx()];

	ncpp::read_data(ncid, cf_latitude, vid_lat, lat);
	ncpp::read_data(ncid, cf_longitude, vid_lon, lon);

	for (size_t ii=0; ii<latlen; ++ii) {
		if (lat_array[ii] != lat[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}

	for (size_t ii=0; ii<lonlen; ++ii) {
		if (lon_array[ii] != lon[ii])
			throw IOException("Error while writing grid - grid and lat/lon coordinates are inconsistent", AT);
	}

	delete[] lat; delete[] lon;
}

void CNRMIO::get_meta_data_ids(const int& ncid, std::map<std::string, int>& map_vid)
{
	const string names[] = {cnrm_altitude, cnrm_latitude, cnrm_longitude, cnrm_aspect, cnrm_slope};
	vector<string> varname(names, names + sizeof(names) / sizeof(names[0]));

	vector<string> dimensions;
	dimensions.push_back(cnrm_points); // All variables have to have the dimension cnrm_points

	for (vector<string>::const_iterator it = varname.begin(); it != varname.end(); ++it) {
		int varid;
		const string& name = *it;

		ncpp::get_variable(ncid, name, varid);
		if (!ncpp::check_dimensions(ncid, name, varid, dimensions))
			throw IOException("Variable '" + name  + "' fails dimension check", AT);

		map_vid[name] = varid;
	}
}

} //namespace
