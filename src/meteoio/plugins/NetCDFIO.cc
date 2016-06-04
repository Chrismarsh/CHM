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
#include "NetCDFIO.h"
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/Timer.h>
#include <meteoio/MathOptim.h>
#include <meteoio/plugins/libncpp.h>

#include <cmath>
#include <cstdio>
#include <algorithm>

using namespace std;

namespace mio {
/**
 * @page netcdf NetCDF
 * @section netcdf_format Format
 * In order to promote creation, access and sharing of scientific data, the NetCDF format has been
 * created as a machine-independent format. NetCDF (network Common Data Form) is therefore an interface
 * for array-oriented data access and a library that provides an implementation of the interface. The
 * <A HREF="http://www.unidata.ucar.edu/downloads/netcdf/index.jsp">NetCDF software</A> was developed
 * at the <A HREF="http://www.unidata.ucar.edu/">Unidata Program Center</A> in Boulder, Colorado.
 * In order to graphicaly explore the content and structure of NetCDF files, you can use the
 * <A HREF="http://www.epic.noaa.gov/java/ncBrowse/">ncBrowse</A> java software.
 *
 * The NetCDF format does not impose a specific set of metadata and therefore in order to easily exchange data
 * within a given field, it is a good idea to standardize the metadata. Several such metadata schema can be used
 * by this plugin:
 * - CF1 - the <A HREF="http://cfconventions.org">conventions</A> for climate and forecast (CF) metadata;
 * - ECMWF - from the <A HREF="http://www.ecmwf.int/">European Centre for Medium-Range Weather Forecasts</A>;
 * - CNRM - from the <A HREF="http://www.cnrm.meteo.fr/">National Centre for Meteorological Research</A>.
 *
 * @section netcdf_units Units
 *
 *
 * @section netcdf_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - DEMFILE: The filename of the file containing the DEM; [Input] section
 * - DEMVAR: The variable name of the DEM within the DEMFILE; [Input] section
 * - GRID2DFILE: the NetCDF file which shall be used for gridded input/output; [Input] and [Output] section
 * - NETCDF_SCHEMA: the schema to use (either CF1 or CNRM or ECMWF); [Input] and [Output] section
 * - DEM_FROM_PRESSURE: if no dem is found but local and sea level pressure grids are found, use them to rebuild a DEM; [Input] section
 *
 * @section netcdf_example Example use
 * @code
 * [Input]
 * DEM     = NETCDF
 * DEMFILE = ./input/Aster_tile.nc
 * @endcode
 *
 * @section netcdf_compilation Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 */

const double NetCDFIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value
const double NetCDFIO::epsilon = 1.0e-10; //when comparing timestamps

const std::string NetCDFIO::cf_time = "time";
const std::string NetCDFIO::cf_latitude = "latitude";
const std::string NetCDFIO::cf_longitude = "longitude";
const std::string NetCDFIO::cf_altitude = "z";

NetCDFIO::NetCDFIO(const std::string& configfile) : cfg(configfile), in_attributes(), out_attributes(), coordin(), coordinparam(), coordout(), coordoutparam(),
                                                    in_dflt_TZ(0.), out_dflt_TZ(0.), in_time_offset(IOUtils::nodata), in_time_multiplier(IOUtils::nodata), 
                                                    dem_altimeter(false), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) : cfg(cfgreader), in_attributes(), out_attributes(), coordin(), coordinparam(), coordout(), coordoutparam(),
                                              in_dflt_TZ(0.), out_dflt_TZ(0.), in_time_offset(IOUtils::nodata), in_time_multiplier(IOUtils::nodata), 
                                              dem_altimeter(false), in_strict(false), out_strict(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::~NetCDFIO() throw() {}

void NetCDFIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("DEM_FROM_PRESSURE", "Input", dem_altimeter, IOUtils::nothrow);
	
	string in_schema = "ECMWF";
	in_schema = IOUtils::strToUpper( cfg.get("NETCDF_SCHEMA", "Input", IOUtils::nothrow) );
	initAttributesMap(in_schema, in_attributes);
	string out_schema = "ECMWF";
	out_schema = IOUtils::strToUpper( cfg.get("NETCDF_SCHEMA", "Output", IOUtils::nothrow) );
	initAttributesMap(out_schema, out_attributes);
}

void NetCDFIO::initAttributesMap(const std::string& schema, std::map<MeteoGrids::Parameters, attributes> &attr)
{
	if (schema.empty()) return;
	
	if (schema=="CF1") {
		attr[MeteoGrids::DEM] = attributes("z", "altitude", "height above mean sea level", "m", IOUtils::nodata);
		attr[MeteoGrids::TA] = attributes("temperature", "air_temperature", "near surface air temperature", "K", IOUtils::nodata);
		attr[MeteoGrids::RH] = attributes("humidity", "relative humidity", "relative humidity", "fraction", IOUtils::nodata);
		attr[MeteoGrids::P] = attributes("pressure", "air_pressure", "near surface air pressure", "Pa", IOUtils::nodata);
	} else if (schema=="CNRM") {
		attr[MeteoGrids::DEM] = attributes("ZS", "", "altitude", "m", IOUtils::nodata);
		attr[MeteoGrids::SLOPE] = attributes("slope", "", "slope angle", "degrees from horizontal", IOUtils::nodata);
		attr[MeteoGrids::AZI] = attributes("aspect", "", "slope aspect", "degrees from north", IOUtils::nodata);
		attr[MeteoGrids::TA] = attributes("Tair", "", "Near Surface Air Temperature", "K", IOUtils::nodata);
		attr[MeteoGrids::RH] = attributes("HUMREL", "", "Relative Humidity", "%", IOUtils::nodata);
		attr[MeteoGrids::QI] = attributes("Qair", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::VW] = attributes("Wind", "", "Wind Speed", "m/s", IOUtils::nodata);
		attr[MeteoGrids::DW] = attributes("Wind_DIR", "", "Wind Direction", "deg", IOUtils::nodata);
		attr[MeteoGrids::PSUM_L] = attributes("Rainf", "", "Rainfall Rate", "kg/m2/s", IOUtils::nodata);
		attr[MeteoGrids::PSUM_S] = attributes("Snowf", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::ISWR_DIR] = attributes("DIR_SWdown", "", "Surface Incident Direct Shortwave Radiation", "W/m2", IOUtils::nodata);
		attr[MeteoGrids::ISWR_DIFF] = attributes("SCA_SWdown", "", "", "", IOUtils::nodata);
		attr[MeteoGrids::P] = attributes("PSurf", "", "Surface Pressure", "Pa", IOUtils::nodata);
		attr[MeteoGrids::ILWR] = attributes("LWdown", "", "Surface Incident Longwave Radiation", "W/m2", IOUtils::nodata);
	} else if (schema=="ECMWF") {
		attr[MeteoGrids::DEM] = attributes("z", "geopotential_height", "geopotential_height", "m", IOUtils::nodata);
		attr[MeteoGrids::TA] = attributes("t2m", "", "2 metre temperature", "K", 2.);
		attr[MeteoGrids::TD] = attributes("d2m", "", "2 metre dewpoint temperature", "K", 2.);
		attr[MeteoGrids::P] = attributes("sp", "surface_air_pressure", "Surface pressure", "Pa", IOUtils::nodata);
		attr[MeteoGrids::P_SEA] = attributes("msl", "air_pressure_at_sea_level", "Mean sea level pressure", "Pa", IOUtils::nodata);
		attr[MeteoGrids::ISWR] = attributes("ssrd", "surface_downwelling_shortwave_flux_in_air", "Surface solar radiation downwards", "J m**-2", IOUtils::nodata);
		attr[MeteoGrids::ILWR] = attributes("strd", "", "Surface thermal radiation downwards", "J m**-2", IOUtils::nodata);
		attr[MeteoGrids::PSUM] = attributes("tp", "", "Total precipitation", "m", IOUtils::nodata);
		attr[MeteoGrids::U] = attributes("u10", "", "10 metre U wind component", "m s**-1", 10.);
		attr[MeteoGrids::V] = attributes("v10", "", "10 metre V wind component", "m s**-1", 10.);
		attr[MeteoGrids::SWE] = attributes("sd", "lwe_thickness_of_surface_snow_amount", "Snow depth", "m of water equivalent", IOUtils::nodata);
		attr[MeteoGrids::TSS] = attributes("skt", "", "Skin temperature", "K", IOUtils::nodata);
	} else
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema+"\"", AT);
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const std::string& arguments)
{
	vector<string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() == 2) {
		read2DGrid_internal(grid_out, vec_argument[0], vec_argument[1]);
	} else {
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::read2DGrid is filename:varname", AT);
	}
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	grid_out.clear();
	const string filename = cfg.get("GRID2DFILE", "Input"); //HACK: also allow using GRID2DPATH
	
	if (read2DGrid_internal(grid_out, filename, parameter, date)) return; //schema naming
	
	if (parameter==MeteoGrids::VW || parameter==MeteoGrids::DW) {	//VW, DW
		Grid2DObject U,V;
		const bool hasU = read2DGrid_internal(U, filename, MeteoGrids::U, date);
		const bool hasV = read2DGrid_internal(V, filename, MeteoGrids::V, date);
		if (hasU==true && hasV==true) {
			grid_out.set(U, IOUtils::nodata);
			if (parameter==MeteoGrids::VW) {
				for (size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
					grid_out(ii) = sqrt( Optim::pow2(U(ii)) + Optim::pow2(V(ii)) );
			} else {
				for (size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
					grid_out(ii) =  fmod( atan2( U(ii), V(ii) ) * Cst::to_deg + 360., 360.); // turn into degrees [0;360)
			}
			return;
		}
	}
	
	if (parameter==MeteoGrids::RH) {								//RH
		Grid2DObject ta;
		DEMObject dem;
		bool hasDEM = false;
		try {
			readDEM(dem);
			hasDEM = true;
		} catch(...){}
		const bool hasQI = read2DGrid_internal(grid_out, filename, MeteoGrids::QI, date);
		const bool hasTA = read2DGrid_internal(ta, filename, MeteoGrids::TA, date);
		if (hasQI && hasDEM && hasTA) {
			for(size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
				grid_out(ii) = Atmosphere::specToRelHumidity(dem(ii), ta(ii), grid_out(ii));
			return;
		}
		
		const bool hasTD = read2DGrid_internal(grid_out, filename, MeteoGrids::TD, date);
		if (hasTA && hasTD) {
			for(size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
				grid_out(ii) = Atmosphere::DewPointtoRh(grid_out(ii), ta(ii), false);
			return;
		}
	}
	
	if (parameter==MeteoGrids::ISWR) {								//ISWR
		Grid2DObject iswr_diff;
		const bool hasISWR_DIFF = read2DGrid_internal(iswr_diff, filename, MeteoGrids::ISWR_DIFF);
		const bool hasISWR_DIR = read2DGrid_internal(grid_out, filename, MeteoGrids::ISWR_DIR);
		if (hasISWR_DIFF && hasISWR_DIR) {
			grid_out += iswr_diff;
			return;
		}
	}
	
	if (parameter==MeteoGrids::PSUM) {								//PSUM
		Grid2DObject psum_s;
		const bool hasPSUM_S = read2DGrid_internal(psum_s, filename, MeteoGrids::PSUM_S);
		const bool hasPSUM_L = read2DGrid_internal(grid_out, filename, MeteoGrids::PSUM_L);
		if (hasPSUM_S && hasPSUM_L) {
			grid_out += psum_s;
			return;
		}
	}
	
	if (parameter==MeteoGrids::PSUM_PH) {								//PSUM_PH
		Grid2DObject psum_l;
		const bool hasPSUM_S = read2DGrid_internal(grid_out, filename, MeteoGrids::PSUM_S);
		const bool hasPSUM_L = read2DGrid_internal(psum_l, filename, MeteoGrids::PSUM_L);
		if (hasPSUM_S && hasPSUM_L) {
			grid_out += psum_l;
			const size_t nrCells = grid_out.getNx()*grid_out.getNy();
			for (size_t ii=0; ii<nrCells; ii++) {
				const double psum = grid_out(ii);
				if (psum!=IOUtils::nodata && psum>0)
					grid_out(ii) = psum_l(ii) / psum;
			}
			return;
		}
	}
	
	throw InvalidArgumentException("Parameter \'"+MeteoGrids::getParameterName(parameter)+"\' either not found in file \'"+filename+"\' or not found in current NetCDF schema", AT);
}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	const string filename = cfg.get("DEMFILE", "Input");
	const string varname = cfg.get("DEMVAR", "Input", IOUtils::nothrow);
	if (!varname.empty()) {
		if (!read2DGrid_internal(dem_out, filename, varname))
			throw InvalidArgumentException("Variable \'"+varname+"\' not found in file \'"+filename+"\'", AT);
	} else {
		if (read2DGrid_internal(dem_out, filename, MeteoGrids::DEM)) return; //schema naming
		if (read2DGrid_internal(dem_out, filename, "Band1")) return; //ASTER naming
		if (read2DGrid_internal(dem_out, filename, "z")) return; //GDAL naming
		
		//last chance: read from pressure grids
		if (dem_altimeter) {
			Grid2DObject p, ta, p_sea;
			if (read2DGrid_internal(p_sea, filename, MeteoGrids::P_SEA) &&
			read2DGrid_internal(p, filename, MeteoGrids::P) && 
			read2DGrid_internal(ta, filename, MeteoGrids::TA)) {
				dem_out.set(p, IOUtils::nodata);
				const double k = Cst::gravity / (Cst::dry_adiabatique_lapse_rate * Cst::gaz_constant_dry_air);
				const double k_inv = 1./k;
				for(size_t ii=0; ii<(dem_out.getNx()*dem_out.getNy()); ii++) {
					const double K = pow(p(ii)/p_sea(ii), k_inv);
					dem_out(ii) = ta(ii)*Cst::earth_R0*(1.-K) / (Cst::dry_adiabatique_lapse_rate * Cst::earth_R0 - ta(ii)*(1.-K));
				}
				return;
			}
		}
		
		throw InvalidArgumentException("The variable containing the DEM could not be found. Please specify it using the DEMVAR key.", AT);
	}
}

void NetCDFIO::readLanduse(Grid2DObject& /*landuse_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readAssimilationData(const Date& /*date_in*/, Grid2DObject& /*da_out*/)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readStationData(const Date&, std::vector<StationData>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readMeteoData(const Date&, const Date&, std::vector< std::vector<MeteoData> >&, const size_t&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::writeMeteoData(const std::vector< std::vector<MeteoData> >&, const std::string&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::readPOI(std::vector<Coords>&)
{
	//Nothing so far
	throw IOException("Nothing implemented here", AT);
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const std::string& arguments)
{
	// arguments is a string of the format filname:varname
	vector<string> vec_argument;
	if (IOUtils::readLineToVec(arguments, vec_argument, ':')  != 2)
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::write2DGrid is filename:varname", AT);

	const string name = vec_argument[1];
	const attributes attr(name, name, name, "", IOUtils::nodata);
	write2DGrid_internal(grid_in, vec_argument[0], attr);
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const string filename = cfg.get("GRID2DFILE", "Output");
	
	const std::map<MeteoGrids::Parameters, attributes>::const_iterator it = in_attributes.find(parameter);
	if (it!=in_attributes.end()) {
		const bool isPrecip = (parameter==MeteoGrids::PSUM || parameter==MeteoGrids::SWE);
		write2DGrid_internal(grid_in, filename, it->second, date, isPrecip);
	} else {
		const string name = MeteoGrids::getParameterName(parameter);
		const attributes attr(name, name, name, "", IOUtils::nodata);
		write2DGrid_internal(grid_in, filename, attr, date);
	}
}

bool NetCDFIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& filename, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const std::map<MeteoGrids::Parameters, attributes>::const_iterator it = in_attributes.find(parameter);
	if (it==in_attributes.end()) return false;
	
	const bool isPrecip = (parameter==MeteoGrids::PSUM || parameter==MeteoGrids::SWE);
	const bool status = read2DGrid_internal(grid_out, filename, it->second.var, date, isPrecip);
	
	return status;
}

//isPrecip is used to convert the units of precip as well as prevent very small precip amounts
bool NetCDFIO::read2DGrid_internal(Grid2DObject& grid_out, const std::string& filename, const std::string& varname, const Date& date, const bool& isPrecip)
{
	int ncid, varid;
	vector<int> dimid, dim_varid;
	vector<string> dimname;
	vector<size_t> dimlen;

	if (!IOUtils::fileExists(filename)) throw FileAccessException(filename, AT); //prevent invalid filenames
	ncpp::open_file(filename, NC_NOWRITE, ncid);
	if (!ncpp::check_variable(ncid, varname)) return false;
	ncpp::get_variable(ncid, varname, varid);
	ncpp::get_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);

	size_t time_index = IOUtils::npos, lat_index = IOUtils::npos, lon_index = IOUtils::npos;
	for(size_t ii=0; ii<dimname.size(); ii++) {
		const string name=dimname[ii];
		if (name=="latitude" || name=="lat")
			lat_index = ii;
		else if (name=="longitude" || name=="lon")
			lon_index = ii;
		else if (name=="time")
			time_index = ii;
		else
			throw InvalidArgumentException("Unknown dimension \'"+name+"\' found in file \'"+filename+"\'", AT);
	}
	if (lat_index==IOUtils::npos || lon_index==IOUtils::npos)
		throw InvalidArgumentException("Both latitudes and longitudes must be provided in file \'"+filename+"\'!", AT);
	if (dimlen[lat_index]<2 || dimlen[lon_index]<2)
		throw IOException("All dimensions for variable '" + varname + "' have to at least have length 2", AT);
	
	//read latitude and longitude vectors
	double *lat = new double[dimlen[lat_index]];
	double *lon = new double[dimlen[lon_index]];
	ncpp::read_data(ncid, dimname[lat_index], dim_varid[lat_index], lat);
	ncpp::read_data(ncid, dimname[lon_index], dim_varid[lon_index], lon);
	
	//read gridded data
	const bool date_requested = (date!=Date());
	const bool date_in_file = (time_index!=IOUtils::npos);
	
	double *grid = new double[dimlen[lat_index]*dimlen[lon_index]];
	if (!date_in_file) {
		if (!date_requested)
			ncpp::read_data(ncid, varname, varid, grid);
		else
			throw InvalidArgumentException("File \'"+filename+"\' does not contain a time dimension", AT);
	} else {
		if (date_requested) {
			const int timeID = dimid[time_index];
			if (in_time_offset==IOUtils::nodata || in_time_multiplier==IOUtils::nodata) {
				getTimeTransform(ncid, timeID, in_time_offset, in_time_multiplier);
			}
			const double timestamp = (date.getJulian() - in_time_offset) / in_time_multiplier;
			const size_t pos = ncpp::find_record(ncid, NetCDFIO::cf_time, timeID, timestamp);
			if (pos == IOUtils::npos) {
				double min, max;
				const bool status = ncpp::get_recordMinMax(ncid, NetCDFIO::cf_time, timeID, min, max);
				if (status) {
					Date d_min, d_max;
					d_min.setDate(min*in_time_multiplier + in_time_offset, in_dflt_TZ);
					d_max.setDate(max*in_time_multiplier + in_time_offset, in_dflt_TZ);
					throw IOException("No record for date " + date.toString(Date::ISO) + ". Records between " + d_min.toString(Date::ISO) + " and " + d_max.toString(Date::ISO), AT);
				}
				else 
					throw IOException("No record for date " + date.toString(Date::ISO), AT);
			}
			ncpp::read_data(ncid, varname, varid, pos, dimlen[lat_index], dimlen[lon_index], grid);
		} else {
			const size_t pos = 1;
			ncpp::read_data(ncid, varname, varid, pos, dimlen[lat_index], dimlen[lon_index], grid);
		}
	}
	
	//read nodata value
	double missing_value=plugin_nodata;
	if (ncpp::check_attribute(ncid, varid, "missing_value")) ncpp::get_attribute(ncid, varname, varid, "missing_value", missing_value);

	//fill our Grid2DObject with all the data that has been read
	ncpp::copy_grid(coordin, coordinparam, dimlen[lat_index], dimlen[lon_index], lat, lon, grid, missing_value, grid_out);
	delete[] lat; delete[] lon; delete[] grid;

	//handle data packing if necessary
	if (ncpp::check_attribute(ncid, varid, "scale_factor")) {
		double scale_factor=1.;
		ncpp::get_attribute(ncid, varname, varid, "scale_factor", scale_factor);
		grid_out *= scale_factor;
	}
	if (ncpp::check_attribute(ncid, varid, "add_offset")) {
		double add_offset=0.;
		ncpp::get_attribute(ncid, varname, varid, "add_offset", add_offset);
		grid_out += add_offset;
	}
	
	//Correct the units if necessary
	if (ncpp::check_attribute(ncid, varid, "units")) {
		string units;
		ncpp::get_attribute(ncid, varname, varid, "units", units);
		if (units=="m**2 s**-2") grid_out /= Cst::gravity;
		if (units=="%") grid_out /= 100.;
		if (units=="J m**-2") grid_out /= (3600.*3.);
		if (isPrecip && units=="m") grid_out *= 100.;
		if (isPrecip) {//reset very low precip to zero
			for (size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
				if (grid_out(ii)<1e-3) grid_out(ii)=0.;
		}
	}

	ncpp::close_file(filename, ncid);
	return true;
}

void NetCDFIO::write2DGrid_internal(Grid2DObject grid_in, const std::string& filename, const attributes& attr, const Date& date, const bool& isPrecip)
{
	const string varname = attr.var;
	const bool is_record = (date != Date() && date!=Date(0.));

	double *lat_array = new double[grid_in.getNy()];
	double *lon_array = new double[grid_in.getNx()];
	int *data = new int[grid_in.getNy() * grid_in.getNx()];
	
	//Correct the units if necessary
	const string units = attr.units;
	if (units=="m**2 s**-2") grid_in *= Cst::gravity;
	if (units=="%") grid_in *= 100.;
	if (units=="J m**-2") grid_in *= (3600.*1.); //HACK: assuming that we do hourly outputs
	if (isPrecip && units=="m") grid_in *= 0.01;
	
	//Compute data packing
	const double data_min = grid_in.grid2D.getMin();
	const double range = grid_in.grid2D.getMax() - data_min;
	double add_offset = 0., scale_factor = 1.;
	if (range!=0.) {
		const long double type_range = 0.5* (static_cast<long double>(std::numeric_limits<int>::max()) - static_cast<long double>(std::numeric_limits<int>::min()));
		scale_factor = static_cast<double> ( range / (type_range - 1.) ); //we reserve the max for nodata
		add_offset = data_min + range/2.; //center the data on the central value of the type
		grid_in -= add_offset;
		grid_in /= scale_factor;
	}
	const double nodata_out = (range!=0)? std::numeric_limits<int>::max() : IOUtils::nodata;
	ncpp::calculate_dimensions(grid_in, lat_array, lon_array);
	ncpp::fill_grid_data(grid_in, nodata_out, data);

	int ncid, did_lat, did_lon, did_time, vid_lat, vid_lon, vid_var, vid_time;
	bool create_dimensions(false), create_variable(false), create_time(false);

	if ( IOUtils::fileExists(filename) ) {
		ncpp::open_file(filename, NC_WRITE, ncid);

		//check of lat/lon are defined and consistent
		if (ncpp::check_dim_var(ncid, cf_latitude) && ncpp::check_dim_var(ncid, cf_longitude)) {
			check_consistency(ncid, grid_in, lat_array, lon_array, did_lat, did_lon, vid_lat, vid_lon);
		} else {
			create_dimensions = true;
		}

		if (is_record) {
			//check if a time dimension/variable already exists
			if (ncpp::check_dim_var(ncid, NetCDFIO::cf_time)) {
				ncpp::get_dimension(ncid, NetCDFIO::cf_time, did_time);
				ncpp::get_variable(ncid, NetCDFIO::cf_time, vid_time);
			} else {
				create_time = true;
			}
		}

		if (ncpp::check_variable(ncid, attr.var)) { // variable exists
			ncpp::get_variable(ncid, attr.var, vid_var);

			vector<int> dimid, dim_varid;
			vector<string> dimname;
			vector<size_t> dimlen;
			ncpp::get_dimension(ncid, attr.var, vid_var, dimid, dim_varid, dimname, dimlen);

			if (is_record) {
				if ((dimname.size() != 3) || (dimname[0] != cf_time) || (dimname[1] != cf_latitude) || (dimname[2] != cf_longitude) || (dimlen[1]!=grid_in.getNy()) || (dimlen[2]!=grid_in.getNx()))
					throw IOException("Variable '" + attr.var  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
			} else {
				if ((dimname[0] != cf_latitude) || (dimname[1] != cf_longitude) || (dimlen[0]!=grid_in.getNy()) || (dimlen[1]!=grid_in.getNx()))
					throw IOException("Variable '" + attr.var  + "' already defined with different dimensions in file '"+ filename  +"'", AT);
			}
		} else {
			create_variable = true;
		}

		ncpp::start_definitions(filename, ncid);
	} else {
		if (!IOUtils::validFileAndPath(filename)) throw InvalidFileNameException(filename, AT);
		ncpp::create_file(filename, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.3");
		create_variable = create_dimensions = true;
		if (is_record) create_time = true;
	}

	if (create_dimensions) create_latlon_dimensions(ncid, grid_in, did_lat, did_lon, vid_lat, vid_lon);
	if (create_time) create_time_dimension(ncid, did_time, vid_time);

	if (is_record && create_variable) {
		ncpp::add_3D_variable(ncid, attr.var, NC_INT, did_time, did_lat, did_lon, vid_var); //NC_DOUBLE or NC_INT or NC_SHORT
		add_attributes_for_variable(ncid, vid_var, attr, nodata_out);
	} else if (create_variable) {
		ncpp::add_2D_variable(ncid, attr.var, NC_INT, did_lat, did_lon, vid_var); //NC_DOUBLE
		add_attributes_for_variable(ncid, vid_var, attr, nodata_out);
	}

	if (range!=0.) {
		ncpp::add_attribute(ncid, vid_var, "add_offset", add_offset);
		ncpp::add_attribute(ncid, vid_var, "scale_factor", scale_factor);
	}
	
	ncpp::end_definitions(filename, ncid);

	if (create_dimensions) {
		ncpp::write_data(ncid, cf_latitude, vid_lat, lat_array);
		ncpp::write_data(ncid, cf_longitude, vid_lon, lon_array);
	}

	if (is_record) {
		size_t pos_start = ncpp::add_record(ncid, NetCDFIO::cf_time, vid_time, static_cast<double>( date.getUnixDate()/3600 ));
		ncpp::write_data(ncid, attr.var, vid_var, grid_in.getNy(), grid_in.getNx(), pos_start, data);
	} else {
		ncpp::write_data(ncid, attr.var, vid_var, data);
	}

	ncpp::close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

void NetCDFIO::getTimeTransform(const int& ncid, const int& varid, double &time_offset, double &time_multiplier) const
{
	string time_units;
	ncpp::get_attribute(ncid, NetCDFIO::cf_time, varid, "units", time_units);
	
	std::vector<std::string> vecString;
	const size_t nrWords = IOUtils::readLineToVec(time_units, vecString);
	if (nrWords<3 || nrWords>4) throw InvalidArgumentException("Invalid format for time units: \'"+time_units+"\'", AT);
	
	if (vecString[0]=="days") time_multiplier = 1.;
	else if (vecString[0]=="hours") time_multiplier = 1./24.;
	else if (vecString[0]=="minutes") time_multiplier = 1./(24.*60.);
	else if (vecString[0]=="seconds") time_multiplier = 1./(24.*3600);
	else throw InvalidArgumentException("Unknown time unit \'"+vecString[0]+"\'", AT);
	
	const string ref_date_str = (nrWords==3)? vecString[2] : vecString[2]+"T"+vecString[3];
	Date refDate;
	if (!IOUtils::convertString(refDate, ref_date_str, in_dflt_TZ))
		throw InvalidArgumentException("Invalid reference date \'"+ref_date_str+"\'", AT);
	
	time_offset = refDate.getJulian();
}

void NetCDFIO::create_latlon_dimensions(const int& ncid, const Grid2DObject& grid_in, int& did_lat, int& did_lon, int& varid_lat, int& varid_lon)
{
	ncpp::add_dimension(ncid, cf_latitude, grid_in.getNy(), did_lat);
	ncpp::add_1D_variable(ncid, cf_latitude, NC_DOUBLE, did_lat, varid_lat);
	ncpp::add_attribute(ncid, varid_lat, "standard_name", cf_latitude);
	ncpp::add_attribute(ncid, varid_lat, "long_name", cf_latitude);
	ncpp::add_attribute(ncid, varid_lat, "units", "degrees_north");

	ncpp::add_dimension(ncid, cf_longitude, grid_in.getNx(), did_lon);
	ncpp::add_1D_variable(ncid, cf_longitude, NC_DOUBLE, did_lon, varid_lon);
	ncpp::add_attribute(ncid, varid_lon, "standard_name", cf_longitude);
	ncpp::add_attribute(ncid, varid_lon, "long_name", cf_longitude);
	ncpp::add_attribute(ncid, varid_lon, "units", "degrees_east");
}

void NetCDFIO::create_time_dimension(const int& ncid, int& did_time, int& varid_time)
{
	ncpp::add_dimension(ncid, NetCDFIO::cf_time, NC_UNLIMITED, did_time);
	ncpp::add_1D_variable(ncid, NetCDFIO::cf_time, NC_DOUBLE, did_time, varid_time);
	ncpp::add_attribute(ncid, varid_time, "standard_name", cf_time);
	ncpp::add_attribute(ncid, varid_time, "long_name", cf_time);
	ncpp::add_attribute(ncid, varid_time, "units", "hours since 1970-1-1");
	ncpp::add_attribute(ncid, varid_time, "calendar", "gregorian");
}

void NetCDFIO::add_attributes_for_variable(const int& ncid, const int& varid, const attributes& attr, const double& nodata_out)
{
	ncpp::add_attribute(ncid, varid, "standard_name", attr.standard_name);
	ncpp::add_attribute(ncid, varid, "long_name", attr.long_name);
	ncpp::add_attribute(ncid, varid, "units", attr.units);
	if (attr.var=="z" || attr.var=="ZS") { //for DEM
		ncpp::add_attribute(ncid, varid, "positive", "up");
		ncpp::add_attribute(ncid, varid, "axis", "Z");
	}
	ncpp::add_attribute(ncid, varid, "missing_value", nodata_out);
}

void NetCDFIO::check_consistency(const int& ncid, const Grid2DObject& grid, double*& lat_array, double*& lon_array,
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

} //namespace
