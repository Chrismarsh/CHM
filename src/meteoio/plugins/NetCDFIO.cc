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
#include <meteoio/plugins/NetCDFIO.h>
#include <meteoio/ResamplingAlgorithms2D.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/FileUtils.h>
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
 * <A HREF="http://www.epic.noaa.gov/java/ncBrowse/">ncBrowse</A> java software or 
 * <A HREF="http://meteora.ucsd.edu/~pierce/ncview_home_page.html">ncview</A>. It is also possible to run *ncdump* on a given
 * file in order to have a look at its structure (such as *ncdump {my_netcdf_file} | more*) and specially the parameters names 
 * (this is useful if remapping is needed, see below for in the \ref netcdf_keywords "keywords" section).
 *
 * The NetCDF format does not impose a specific set of metadata and therefore in order to easily exchange data
 * within a given field, it is a good idea to standardize the metadata. Several such metadata schema can be used
 * by this plugin:
 * - CF1 - the <A HREF="http://cfconventions.org">conventions</A> for climate and forecast (CF) metadata;
 * - ECMWF - from the <A HREF="http://www.ecmwf.int/">European Centre for Medium-Range Weather Forecasts</A>, see the <A HREF="https://software.ecmwf.int/wiki/display/TIGGE/Soil+temperature">ECMWF Wiki</A> for a description of the available fields;
 * - CNRM - from the <A HREF="http://www.cnrm.meteo.fr/">National Centre for Meteorological Research</A>;
 * - WRF - the <A HREF="http://www.wrf-model.org/index.php">Weather Research & Forecasting</A> model.
 * 
 * @section netcdf_compilation Compilation
 * In order to compile this plugin, you need libnetcdf (for C). For Linux, please select both the libraries and
 * their development files in your package manager.
 *
 * @section netcdf_keywords Keywords
 * This plugin uses the following keywords:
 * - COORDSYS: coordinate system (see Coords); [Input] and [Output] section
 * - COORDPARAM: extra coordinates parameters (see Coords); [Input] and [Output] section
 * - DEMFILE: The filename of the file containing the DEM; [Input] section
 * - DEMVAR: The variable name of the DEM within the DEMFILE; [Input] section
 * - GRID2DPATH: if this directory contains files, they will be used for reading the input from; [Input] section
 * - METEO_EXT: only the files containing this pattern in their filename will be used; [Input] section
 * - GRID2DFILE: if GRID2DPATH has not been defined or if it does not contain files matching the METEO_EXT extension, provides
 * the NetCDF file which shall be used for gridded input/output; [Input] and [Output] section
 * - NETCDF_SCHEMA: the schema to use (either CF1 or CNRM or ECMWF or WRF); [Input] and [Output] section
 * - NETCDF::{MeteoGrids::Parameters} = {netcdf_param_name} : this allows to remap the names as found in the NetCDF file to the MeteoIO grid parameters; [Input] section;
 * - DEM_FROM_PRESSURE: if no dem is found but local and sea level pressure grids are found, use them to rebuild a DEM; [Input] section
 *
 * When providing multiple files in one directory, in case of overlapping files, the file containing the newest data has priority. This is
 * convenient when using forecats data to automatically use the most short-term forecast.
 * 
 * @section netcdf_example Example use
 * Using this plugin to build downscaled time series at virtual stations, with the ECMWF Era Interim data set (see section below):
 * @code
 * [Input]
 * GRID2D    = NETCDF
 * GRID2DPATH =  /data/meteo_reanalysis
 * METEO_EXT = .nc
 * NETCDF_SCHEMA = ECMWF
 * 
 * DEM = NETCDF
 * DEMFILE = /data/meteo_reanalysis/ECMWF_Europe_20150101-20150701.nc
 * DEM_FROM_PRESSURE = true
 * 
 * #The lines below have nothing to do with this plugin
 * Downscaling = true
 * VSTATION1 = 46.793029 9.821343 ;this is Davos
 * Virtual_parameters = TA RH PSUM ISWR ILWR P VW DW TSS HS RSWR TSG ;this has to fit the parameter set in the data files
 * @endcode
 * 
 * Another example, to extract precipitation from the MeteoSwiss daily precipitation reanalysis, RhiresD
 * @code
 * [Input]
 * DEM     = NETCDF
 * DEMFILE = ./input/ch02_lonlat.nc
 * 
 * GRID2D    = NETCDF
 * GRID2DPATH =  /data/meteo_reanalysis
 * METEO_EXT = .nc
 * NETCDF::PSUM = RhiresD               ;overwrite the PSUM parameter with "RhiresD", for example for MeteoCH reanalysis
 * 
 * #The lines below have nothing to do with this plugin
 * Downscaling = true
 * VSTATION1 = 46.793029 9.821343 ;this is Davos
 * Virtual_parameters = PSUM ;this has to fit the parameter set in the data files
 * @endcode
 * 
 * @section netcdf_meteoch MeteoCH RhiresD & similar products
 * <A HREF="http://www.meteoswiss.admin.ch/home.html?tab=overview">MeteoSwiss</A> provides <A HREF="http://www.ifu.ethz.ch/hydrologie/research/research_data/proddoc.pdf">reanalysis</A> of precipitation and other meteo fields from 1961 to present over Switzerland for different time scales: daily, monthly, yearly, as well as hourly (CPC dataset). The DEM are also provided, either in lat/lon, 
 * Swiss coordinates, rotated lat/lon, ... These data sets must be requested from MeteoSwiss and are available with a specific license for research. 
 * 
 * @section netcdf_wrf WRF output files
 * While <A HREF="http://www.wrf-model.org/index.php">WRF</A> can write its <A HREF="http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap5.htm#fields">outputs in NetCDF</A>, unfortunately
 * it does not follow the CF1 convention and relies on lots of idiosyncracies (see http://www.ncl.ucar.edu/Applications/wrfnetcdf.shtml) that break lots of
 * applications dealing with NetCDF. If some fields are not read by MeteoIO,
 * please follow the tips given \ref netcdf_tricks "below". Moreover, WRF assumes that latitudes / longitudes are given on an ideal sphere while standard
 * coordinates systems assume an ellipsoid. This may lead to trouble when converting model coordinates to real world coordinates (see
 * http://www.pkrc.net/wrf-lambert.html).
 *
 * @section netcdf_ecmwf ECMWF Era Interim
 * The Era Interim data can be downloaded on the <A HREF="http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/">ECMWF dataserver</A> 
 * after creating an account and login in. 
 * 
 * It is recommended to extract data at 00:00, and 12:00 for all steps 3, 6, 9, 12. The select the following fields:
 * 10 metre U wind component, 10 metre V wind component, 2 metre dewpoint temperature, 2 metre temperature, Forecast albedo, Mean sea level pressure, Skin temperature, Snow density, Snow depth, Soil temperature level 1, Surface pressure, Surface solar radiation downwards, Surface thermal radiation downwards, Total precipitation
 * 
 * Here we have included the *forecast albedo* so the RSWR can be computed from ISWR and the *mean sea level pressure* and *surface pressure*
 * as proxies to compute the elevation. If you have the altitude in a separate file, it can be declared as DEM and there would be no need for the sea 
 *level pressure (this would also be much more precise).
 * 
 * You should therefore have the following request:
 * @code
 * Parameter: 10 metre U wind component, 10 metre V wind component, 2 metre dewpoint temperature, 2 metre temperature, Forecast albedo, 
 *            Mean sea level pressure, Skin temperature, Snow density, Snow depth, Soil temperature level 1, Surface pressure, 
 *            Surface solar radiation downwards, Surface thermal radiation downwards, Total precipitation
 *      Step: 3 to 12 by 3
 *      Type: Forecast
 *      Time: 00:00:00, 12:00:00
 * @endcode
 * 
 * With the <A HREF="https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets">ECMWF Python Library</A>, the request 
 * would be for example:
 * @code
 * #!/usr/bin/env python
 * from ecmwfapi import ECMWFDataServer
 * server = ECMWFDataServer()
 * server.retrieve({
 * "class": "ei",
 * "dataset": "interim",
 * "date": "2015-01-01/to/2015-01-31",
 * "expver": "1",
 * "grid": "0.75/0.75",
 * "levtype": "sfc",
 * "param": "33.128/134.128/139.128/141.128/151.128/165.128/166.128/167.128/168.128/169.128/175.128/205.128/228.128/235.128/243.128",
 * "step": "3/6/9/12",
 * "area":"42.2/-1.5/51.7/15.7",
 * "stream": "oper",
 * "format":"netcdf",
 * "target": "my-era-interim.nc",
 * "time": "00/12",
 * "type": "fc",
 * })
 * @endcode
 *
 * @section netcdf_tricks Saving the day when a file is not standard compliant
 * Unfortunatelly, the naming of the parameters and dimensions within the files is not always standard nor consistent. In order to handle the parameters names,
 * simply run *ncdump {my_netcdf_file} | more* and use the name mapping facility of this plugin to map the non-standard parameters to our internal names
 * (see the \ref netcdf_keywords "plugin keywords"). When the dimensions are not standard (for example the time axis being called "TIME_T"),
 * use first the <A HREF="http://linux.die.net/man/1/ncrename">ncrename</A> tool that is part of the
 * <A HREF="http://nco.sourceforge.net/">NCO utilities</A> to rename both the dimension (-d) and the variable (-v):
 * @code
 * ncrename -d TIME_T,time -v TIME_T,time {my_netcdf_file}
 * @endcode
 */
 
const double NetCDFIO::plugin_nodata = -9999999.; //CNRM-GAME nodata value
const std::string NetCDFIO::cf_time = "time";
const std::string NetCDFIO::cf_latitude = "latitude";
const std::string NetCDFIO::cf_longitude = "longitude";

NetCDFIO::NetCDFIO(const std::string& configfile) : cfg(configfile), cache_meteo_files(), in_attributes(), out_attributes(), 
                                                    coordin(), coordinparam(), coordout(), coordoutparam(), time_dimension(),
                                                    in_dflt_TZ(0.), out_dflt_TZ(0.), in_time_offset(0.), in_time_multiplier(1.),
                                                    dem_altimeter(false), in_strict(false), out_strict(false), meteo_cache_ready(false), wrf_hacks(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

NetCDFIO::NetCDFIO(const Config& cfgreader) : cfg(cfgreader), cache_meteo_files(), in_attributes(), out_attributes(), 
                                              coordin(), coordinparam(), coordout(), coordoutparam(), time_dimension(),
                                              in_dflt_TZ(0.), out_dflt_TZ(0.), in_time_offset(0.), in_time_multiplier(1.),
                                              dem_altimeter(false), in_strict(false), out_strict(false), meteo_cache_ready(false), wrf_hacks(false), vecMetaData()
{
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);
	parseInputOutputSection();
}

void NetCDFIO::parseInputOutputSection()
{
	//default timezones
	in_dflt_TZ = out_dflt_TZ = IOUtils::nodata;
	cfg.getValue("TIME_ZONE", "Input", in_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("TIME_ZONE", "Output", out_dflt_TZ, IOUtils::nothrow);
	cfg.getValue("DEM_FROM_PRESSURE", "Input", dem_altimeter, IOUtils::nothrow);
	
	const std::string in_schema( IOUtils::strToUpper( cfg.get("NETCDF_SCHEMA", "Input", IOUtils::nothrow) ) );
	if (!in_schema.empty()) initAttributesMap(in_schema, in_attributes);
	else initAttributesMap("ECMWF", in_attributes);
	
	const std::string out_schema( IOUtils::strToUpper( cfg.get("NETCDF_SCHEMA", "Output", IOUtils::nothrow) ) );
	if (!out_schema.empty()) initAttributesMap(out_schema, out_attributes);
	else initAttributesMap("ECMWF", out_attributes);
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
		attr[MeteoGrids::TSG] = attributes("stl1", "surface_temperature", "Soil temperature level 1", "K", IOUtils::nodata); //this is from 0 to -7cm
		attr[MeteoGrids::ALB] = attributes("al", "surface_albedo", "Albedo", "(0 - 1)", IOUtils::nodata);
		attr[MeteoGrids::ALB] = attributes("fal", "", "Forecast albedo", "(0 - 1)", IOUtils::nodata);
		attr[MeteoGrids::RSNO] = attributes("rsn", "", "Snow density", "kg m**-3", IOUtils::nodata);
		attr[MeteoGrids::ROT] = attributes("ro", "", "Runoff", "m", IOUtils::nodata);
	} else if (schema=="WRF") {
		attr[MeteoGrids::DEM] = attributes("HGT", "Terrain Height", "Terrain Height", "m", IOUtils::nodata);
		attr[MeteoGrids::P] = attributes("PSFC", "Surface pressure", "Surface pressure", "Pa", IOUtils::nodata);
		attr[MeteoGrids::TA] = attributes("T2", "2-meter temperature", "2-meter temperature", "K", 2.);
		attr[MeteoGrids::QI] = attributes("Q2", "2-meter specific humidity", "2-meter specific humidity", "kg kg-1", 2);
		attr[MeteoGrids::ISWR] = attributes("ACSWDNB", "Downward SW surface radiation", "Downward SW surface radiation", "W m**-2", IOUtils::nodata);
		attr[MeteoGrids::RSWR] = attributes("ACSWUPB", "Upwelling Surface Shortwave Radiation", "Upwelling Surface Shortwave Radiation", "W m**-2", IOUtils::nodata);
		attr[MeteoGrids::ILWR] = attributes("ACLWDNB", "Downward LW surface radiation", "Downward LW surface radiation", "W m**-2", IOUtils::nodata);
		attr[MeteoGrids::ROT] = attributes("SFROFF", "Surface runoff ", "Surface runoff ", "kg*m2*s-1", IOUtils::nodata);
		attr[MeteoGrids::HS] = attributes("SNOWH", "Snow depth", "Snow depth", "Pa", IOUtils::nodata);
		attr[MeteoGrids::TSS] = attributes("TSK", "Surface skin temperature", "Surface skin temperature", "K", IOUtils::nodata);
		attr[MeteoGrids::U] = attributes("U10", "10-meter wind speed", "10 metre U wind component", "m s**-1", 10.);
		attr[MeteoGrids::V] = attributes("V10", "10-meter wind speed", "10 metre V wind component", "m s**-1", 10.);
		wrf_hacks = true;
	} else
		throw InvalidArgumentException("Invalid schema selected for NetCDF: \""+schema+"\"", AT);
	
	std::vector<std::string> custom_attr;
	const size_t nrOfCustoms = cfg.findKeys(custom_attr, "NETCDF::", "Input");
	for (size_t ii=0; ii<nrOfCustoms; ++ii) {
		const size_t found = custom_attr[ii].find_last_of(":");
		if (found==std::string::npos || found==custom_attr[ii].length()) continue;

		const std::string meteo_grid( custom_attr[ii].substr(found+1) );
		const std::string netcdf_param = cfg.get(custom_attr[ii], "Input");
		const size_t param_index = MeteoGrids::getParameterIndex(meteo_grid);
		if (param_index==IOUtils::npos)
			throw InvalidArgumentException("Parameter '"+meteo_grid+"' is not a valid MeteoGrid! Please correct key '"+custom_attr[ii]+"'", AT);
		
		attr[ static_cast<MeteoGrids::Parameters>(param_index) ] = attributes(netcdf_param, "", "", "", IOUtils::nodata);
	}
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const std::string& arguments)
{
	std::vector<std::string> vec_argument;
	IOUtils::readLineToVec(arguments, vec_argument, ':');

	if (vec_argument.size() == 2) {
		read2DGrid_internal(grid_out, vec_argument[0], vec_argument[1]);
	} else {
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::read2DGrid is filename:varname", AT);
	}
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (!meteo_cache_ready) {
		const std::string in_grid2d_path = cfg.get("GRID2DPATH", "Input", IOUtils::nothrow);
		if (!in_grid2d_path.empty()) scanMeteoPath(in_grid2d_path,  cache_meteo_files);
		meteo_cache_ready = true;
	}

	if (!cache_meteo_files.empty()) {
		for (size_t ii=0; ii<cache_meteo_files.size(); ii++) {
			const Date date_start( cache_meteo_files[ii].first.first );
			const Date date_end( cache_meteo_files[ii].first.second );
			if (date>=date_start && date<=date_end) {
				const std::string filename( cache_meteo_files[ii].second );
				read2DGrid(grid_out, parameter, date, filename);
				return;
			}
		}
		//the date was not found
		std::string in_grid2d_path;
		cfg.getValue("GRID2DPATH", "Input", in_grid2d_path);
		throw InvalidArgumentException("No Gridded data found for "+date.toString(Date::ISO)+"in '"+in_grid2d_path+"'", AT);
	} else {
		const std::string filename = cfg.get("GRID2DFILE", "Input");
		read2DGrid(grid_out, parameter, date, filename);
	}
}

void NetCDFIO::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date, const std::string& filename)
{
	grid_out.clear();
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
			for (size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
				grid_out(ii) = Atmosphere::specToRelHumidity(dem(ii), ta(ii), grid_out(ii));
			return;
		}
		
		const bool hasTD = read2DGrid_internal(grid_out, filename, MeteoGrids::TD, date);
		if (hasTA && hasTD) {
			for (size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
				grid_out(ii) = Atmosphere::DewPointtoRh(grid_out(ii), ta(ii), false);
			return;
		}
	}
	
	if (parameter==MeteoGrids::RSWR) {								//RSWR
		bool hasISWR = false;
		try {
			read2DGrid(grid_out, MeteoGrids::ISWR, date);
			hasISWR = true;
		} catch(...){}
		Grid2DObject alb;
		const bool hasALB = read2DGrid_internal(alb, filename, MeteoGrids::ALB);
		if (hasALB && hasISWR) {
			grid_out *= alb;
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
	
	if (parameter==MeteoGrids::HS) {								//HS
		Grid2DObject rsno;
		const bool hasRSNO = read2DGrid_internal(rsno, filename, MeteoGrids::RSNO);
		const bool hasSWE = read2DGrid_internal(grid_out, filename, MeteoGrids::SWE);
		if (hasRSNO && hasSWE) {
			grid_out *= 1000.0; //convert mm=kg/m^3 into kg
			grid_out /= rsno;
			return;
		}
	}
	
	throw InvalidArgumentException("Parameter \'"+MeteoGrids::getParameterName(parameter)+"\' either not found in file \'"+filename+"\' or not found in current NetCDF schema", AT);
}

void NetCDFIO::readDEM(DEMObject& dem_out)
{
	const std::string filename = cfg.get("DEMFILE", "Input");
	const std::string varname = cfg.get("DEMVAR", "Input", IOUtils::nothrow);
	if (!varname.empty()) {
		if (!read2DGrid_internal(dem_out, filename, varname))
			throw InvalidArgumentException("Variable \'"+varname+"\' not found in file \'"+filename+"\'", AT);
	} else {
		if (read2DGrid_internal(dem_out, filename, MeteoGrids::DEM)) return; //schema naming
		if (read2DGrid_internal(dem_out, filename, "Band1")) return; //ASTER naming
		if (read2DGrid_internal(dem_out, filename, "z")) return; //GDAL naming
		if (read2DGrid_internal(dem_out, filename, "height")) return; //MeteoCH naming
		if (read2DGrid_internal(dem_out, filename, "HGT")) return; //WRF naming
		
		//last chance: read from pressure grids
		if (dem_altimeter) {
			Grid2DObject p, ta, p_sea;
			if (read2DGrid_internal(p_sea, filename, MeteoGrids::P_SEA) &&
			read2DGrid_internal(p, filename, MeteoGrids::P) && 
			read2DGrid_internal(ta, filename, MeteoGrids::TA)) {
				dem_out.set(p, IOUtils::nodata);
				const double k = Cst::gravity / (Cst::mean_adiabatique_lapse_rate * Cst::gaz_constant_dry_air);
				const double k_inv = 1./k;
				for (size_t ii=0; ii<(dem_out.getNx()*dem_out.getNy()); ii++) {
					const double K = pow(p(ii)/p_sea(ii), k_inv);
					dem_out(ii) = ta(ii)*Cst::earth_R0*(1.-K) / (Cst::mean_adiabatique_lapse_rate * Cst::earth_R0 - ta(ii)*(1.-K));
				}
				return;
			}
		}
		
		throw InvalidArgumentException("The variable containing the DEM could not be found. Please specify it using the DEMVAR key.", AT);
	}
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const std::string& arguments)
{
	// arguments is a string of the format filname:varname
	std::vector<std::string> vec_argument;
	if (IOUtils::readLineToVec(arguments, vec_argument, ':')  != 2)
		throw InvalidArgumentException("The format for the arguments to NetCDFIO::write2DGrid is filename:varname", AT);

	const std::string name( vec_argument[1] );
	const attributes attr(name, name, name, "", IOUtils::nodata);
	write2DGrid_internal(grid_in, vec_argument[0], attr);
}

void NetCDFIO::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	const std::string filename = cfg.get("GRID2DFILE", "Output");
	
	const std::map<MeteoGrids::Parameters, attributes>::const_iterator it = in_attributes.find(parameter);
	if (it!=in_attributes.end()) {
		const bool isPrecip = (parameter==MeteoGrids::PSUM || parameter==MeteoGrids::SWE);
		write2DGrid_internal(grid_in, filename, it->second, date, isPrecip);
	} else {
		const std::string name( MeteoGrids::getParameterName(parameter) );
		const attributes attr(name, name, name, "", IOUtils::nodata);
		write2DGrid_internal(grid_in, filename, attr, date);
	}
}

//custom function for sorting cache_meteo_files
struct sort_pred {
	bool operator()(const std::pair<std::pair<Date,Date>,string> &left, const std::pair<std::pair<Date,Date>,string> &right) {
		if (left.first.first < right.first.first) return true;
		if (left.first.first > right.first.first) return false;
		return left.first.second < right.first.second; //date_start equallity case
	}
};

void NetCDFIO::scanMeteoPath(const std::string& meteopath_in,  std::vector< std::pair<std::pair<mio::Date, mio::Date>, std::string> > &meteo_files)
{
	meteo_files.clear();

	const std::string meteo_ext = cfg.get("METEO_EXT", "INPUT", IOUtils::nothrow);
	std::list<std::string> dirlist = FileUtils::readDirectory(meteopath_in, meteo_ext);
	if (dirlist.empty()) return; //nothing to do if the directory is empty, we will transparently swap to using GRID2DFILE
	dirlist.sort();

	//Check date range in every filename and cache it
	std::list<std::string>::const_iterator it = dirlist.begin();
	while ((it != dirlist.end())) {
		const std::string filename( meteopath_in + "/" + *it );
		
		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		int ncid;
		ncpp::open_file(filename, NC_NOWRITE, ncid);
		if (time_dimension.empty()) ncpp::get_unlimited_dimname(ncid, time_dimension);
		if (!wrf_hacks) getTimeTransform(ncid, in_time_offset, in_time_multiplier); //always re-read offset and multiplier, it might be different for each file
		double min, max;
		const bool status = (!wrf_hacks)? ncpp::get_dimensionMinMax(ncid, time_dimension, min, max) : ncpp::get_wrf_dimensionMinMax(ncid, time_dimension, min, max);
		ncpp::close_file(filename, ncid); //no need to keep file open anymore
		
		if (!status) throw IOException("Could not get min/max time for file '"+filename+"'", AT);
		
		const Date d_min(min*in_time_multiplier + in_time_offset, in_dflt_TZ);
		const Date d_max(max*in_time_multiplier + in_time_offset, in_dflt_TZ);
		const std::pair<Date, Date> dateRange(d_min, d_max);
		const std::pair<std::pair<Date, Date>,std::string> tmp(dateRange, filename);
		meteo_files.push_back(tmp);
		it++;
	}
	std::sort(meteo_files.begin(), meteo_files.end(), sort_pred());
	
	//now handle overlaping files: truncate the end date of the file starting earlier
	for (size_t ii=0; ii<(meteo_files.size()-1); ii++) {
		if (meteo_files[ii].first.second > meteo_files[ii+1].first.first)
			meteo_files[ii].first.second = meteo_files[ii+1].first.first;
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
	if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames

	int ncid, varid;
	ncpp::open_file(filename, NC_NOWRITE, ncid);
	if (!ncpp::check_variable(ncid, varname)) {
		ncpp::close_file(filename, ncid);
		return false;
	}
	ncpp::get_variable(ncid, varname, varid);

	std::vector<int> dimid, dim_varid;
	std::vector<string> dimname;
	std::vector<size_t> dimlen;
	if (!wrf_hacks)
		ncpp::get_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);
	else
		ncpp::get_wrf_dimension(ncid, varname, varid, dimid, dim_varid, dimname, dimlen);

	if (time_dimension.empty()) ncpp::get_unlimited_dimname(ncid, time_dimension);

	size_t time_index = IOUtils::npos, lat_index = IOUtils::npos, lon_index = IOUtils::npos;
	for (size_t ii=0; ii<dimname.size(); ii++) {
		const std::string name( IOUtils::strToLower(dimname[ii]) );
		if (name=="latitude" || name=="lat" || name=="south_north")
			lat_index = ii;
		else if (name=="longitude" || name=="lon" || name=="west_east")
			lon_index = ii;
		else if (name==IOUtils::strToLower(time_dimension))
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
	if (!wrf_hacks) {
		ncpp::read_data(ncid, dimname[lat_index], dim_varid[lat_index], lat);
		ncpp::read_data(ncid, dimname[lon_index], dim_varid[lon_index], lon);
	} else
		ncpp::read_wrf_latlon(ncid, dim_varid[lat_index], dim_varid[lon_index], dimlen[lat_index], dimlen[lon_index], lat, lon);

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
			if (!wrf_hacks) getTimeTransform(ncid, in_time_offset, in_time_multiplier);
			const double timestamp = (date.getJulian() - in_time_offset) / in_time_multiplier;
			const size_t pos = (!wrf_hacks)? ncpp::find_record(ncid, time_dimension, timestamp) : ncpp::find_wrf_record(ncid, time_dimension, timestamp);
			if (pos == IOUtils::npos) {
				double min, max;
				const bool status = (!wrf_hacks)? ncpp::get_dimensionMinMax(ncid, time_dimension, min, max) : ncpp::get_wrf_dimensionMinMax(ncid, time_dimension, min, max);
				if (status) {
					Date d_min, d_max;
					d_min.setDate(min*in_time_multiplier + in_time_offset, in_dflt_TZ);
					d_max.setDate(max*in_time_multiplier + in_time_offset, in_dflt_TZ);
					throw IOException("No record for date " + date.toString(Date::ISO) + ". Records between " + d_min.toString(Date::ISO) + " and " + d_max.toString(Date::ISO)+" in file '"+filename+"'", AT);
				}
				else 
					throw IOException("No record for date " + date.toString(Date::ISO), AT);
			}
			ncpp::read_data(ncid, varname, varid, pos, dimlen[lat_index], dimlen[lon_index], grid);
		} else {
			const size_t pos = 0;
			ncpp::read_data(ncid, varname, varid, pos, dimlen[lat_index], dimlen[lon_index], grid);
		}
	}
	
	//read nodata value
	double missing_value = plugin_nodata;
	if (ncpp::check_attribute(ncid, varid, "missing_value")) ncpp::get_attribute(ncid, varname, varid, "missing_value", missing_value);

	//fill our Grid2DObject with all the data that has been read
	ncpp::copy_grid(coordin, coordinparam, dimlen[lat_index], dimlen[lon_index], lat, lon, grid, missing_value, grid_out);
	delete[] lat; delete[] lon; delete[] grid;

	//handle data packing if necessary
	if (ncpp::check_attribute(ncid, varid, "scale_factor")) {
		double scale_factor = 1.;
		ncpp::get_attribute(ncid, varname, varid, "scale_factor", scale_factor);
		grid_out *= scale_factor;
	}
	if (ncpp::check_attribute(ncid, varid, "add_offset")) {
		double add_offset = 0.;
		ncpp::get_attribute(ncid, varname, varid, "add_offset", add_offset);
		grid_out += add_offset;
	}
	
	//Correct the units if necessary
	if (ncpp::check_attribute(ncid, varid, "units")) {
		std::string units;
		ncpp::get_attribute(ncid, varname, varid, "units", units);
		if (units=="m**2 s**-2") grid_out /= Cst::gravity;
		if (units=="%") grid_out /= 100.;
		if (units=="J m**-2") grid_out /= (3600.*3.);
		if (isPrecip && units=="m") grid_out *= 100.;
		if (isPrecip) {//reset very low precip to zero
			for (size_t ii=0; ii<(grid_out.getNx()*grid_out.getNy()); ii++)
				if (grid_out(ii)<1e-3 && grid_out(ii)!=mio::IOUtils::nodata) grid_out(ii)=0.;
		}
	}

	ncpp::close_file(filename, ncid);
	return true;
}

void NetCDFIO::write2DGrid_internal(Grid2DObject grid_in, const std::string& filename, const attributes& attr, const Date& date, const bool& isPrecip)
{
	const std::string varname( attr.var );
	const bool is_record = (date != Date() && date!=Date(0.));

	double *lat_array = new double[grid_in.getNy()];
	double *lon_array = new double[grid_in.getNx()];
	int *data = new int[grid_in.getNy() * grid_in.getNx()];
	
	//Correct the units if necessary
	const std::string units( attr.units );
	if (units=="m**2 s**-2") grid_in *= Cst::gravity;
	if (units=="%") grid_in *= 100.;
	if (units=="J m**-2") grid_in *= (3600.*1.); //HACK: assuming that we do hourly outputs
	if (isPrecip && units=="m") grid_in *= 0.01;
	
	ncpp::calculate_dimensions(grid_in, lat_array, lon_array);
	ncpp::fill_grid_data(grid_in, IOUtils::nodata, data);

	int ncid, did_lat, did_lon, did_time, vid_lat, vid_lon, vid_var, vid_time;
	bool create_dimensions(false), create_variable(false), create_time(false);

	if ( FileUtils::fileExists(filename) ) {
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

			std::vector<int> dimid, dim_varid;
			std::vector<string> dimname;
			std::vector<size_t> dimlen;
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
		if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException(filename, AT);
		ncpp::create_file(filename, NC_CLASSIC_MODEL, ncid);
		ncpp::add_attribute(ncid, NC_GLOBAL, "Conventions", "CF-1.3");
		create_variable = create_dimensions = true;
		if (is_record) create_time = true;
	}

	if (create_dimensions) create_latlon_dimensions(ncid, grid_in, did_lat, did_lon, vid_lat, vid_lon);
	if (create_time) create_time_dimension(ncid, did_time, vid_time);

	if (is_record && create_variable) {
		ncpp::add_3D_variable(ncid, attr.var, NC_INT, did_time, did_lat, did_lon, vid_var); //NC_DOUBLE or NC_INT or NC_SHORT
		add_attributes_for_variable(ncid, vid_var, attr, IOUtils::nodata);
	} else if (create_variable) {
		ncpp::add_2D_variable(ncid, attr.var, NC_INT, did_lat, did_lon, vid_var); //NC_DOUBLE
		add_attributes_for_variable(ncid, vid_var, attr, IOUtils::nodata);
	}
	ncpp::end_definitions(filename, ncid);

	if (create_dimensions) {
		ncpp::write_data(ncid, cf_latitude, vid_lat, lat_array);
		ncpp::write_data(ncid, cf_longitude, vid_lon, lon_array);
	}

	if (is_record) {
		const size_t pos_start = ncpp::add_record(ncid, NetCDFIO::cf_time, vid_time, static_cast<double>(date.getUnixDate())/3600. );
		ncpp::write_data(ncid, attr.var, vid_var, grid_in.getNy(), grid_in.getNx(), pos_start, data);
	} else {
		ncpp::write_data(ncid, attr.var, vid_var, data);
	}

	ncpp::close_file(filename, ncid);
	delete[] lat_array; delete[] lon_array; delete[] data;
}

void NetCDFIO::getTimeTransform(const int& ncid, double &time_offset, double &time_multiplier) const
{
	if (time_dimension.empty()) throw IOException("time_dimension has not been initialized!", AT);
	const std::string time_units = ncpp::get_DimAttribute(ncid, time_dimension, "units");

	std::vector<std::string> vecString;
	const size_t nrWords = IOUtils::readLineToVec(time_units, vecString);
	if (nrWords<3 || nrWords>4) throw InvalidArgumentException("Invalid format for time units: \'"+time_units+"\'", AT);
	
	const double equinox_year = 365.242198781; //definition used by the NetCDF Udunits package
	
	if (vecString[0]=="years") time_multiplier = equinox_year;
	else if (vecString[0]=="months") time_multiplier = equinox_year/12.;
	else if (vecString[0]=="days") time_multiplier = 1.;
	else if (vecString[0]=="hours") time_multiplier = 1./24.;
	else if (vecString[0]=="minutes") time_multiplier = 1./(24.*60.);
	else if (vecString[0]=="seconds") time_multiplier = 1./(24.*3600);
	else throw InvalidArgumentException("Unknown time unit \'"+vecString[0]+"\'", AT);
	
	const std::string ref_date_str = (nrWords==3)? vecString[2] : vecString[2]+"T"+vecString[3];
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
