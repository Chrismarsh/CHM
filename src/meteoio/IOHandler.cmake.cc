/***********************************************************************************/
/*  Copyright 2009-2012 WSL Institute for Snow and Avalanche Research  SLF-DAVOS   */
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
#include <meteoio/IOHandler.h>
#include <meteoio/IOExceptions.h>
#include <meteoio/IOUtils.h>
#include <meteoio/MathOptim.h>
#include <meteoio/dataClasses/MeteoData.h> //needed for the merge strategies

#include <algorithm>
#include <fstream>

#cmakedefine PLUGIN_ALPUG
#cmakedefine PLUGIN_ARCIO
#cmakedefine PLUGIN_A3DIO
#cmakedefine PLUGIN_ARPSIO
#cmakedefine PLUGIN_CNRMIO
#cmakedefine PLUGIN_DBO
#cmakedefine PLUGIN_GRASSIO
#cmakedefine PLUGIN_GEOTOPIO
#cmakedefine PLUGIN_SMETIO
#cmakedefine PLUGIN_SNIO
#cmakedefine PLUGIN_PGMIO
#cmakedefine PLUGIN_IMISIO
#cmakedefine PLUGIN_OSHDIO
#cmakedefine PLUGIN_GRIBIO
#cmakedefine PLUGIN_PNGIO
#cmakedefine PLUGIN_BORMAIO
#cmakedefine PLUGIN_COSMOXMLIO
#cmakedefine PLUGIN_GSNIO
#cmakedefine PLUGIN_NETCDFIO
#cmakedefine PLUGIN_PSQLIO
#cmakedefine PLUGIN_SASEIO

#include <meteoio/plugins/ALPUG.h>
#include <meteoio/plugins/ARCIO.h>
#include <meteoio/plugins/A3DIO.h>
#include <meteoio/plugins/ARPSIO.h>
#include <meteoio/plugins/GrassIO.h>
#include <meteoio/plugins/GeotopIO.h>
#include <meteoio/plugins/PGMIO.h>
#include <meteoio/plugins/SMETIO.h>
#include <meteoio/plugins/SNIO.h>

#ifdef PLUGIN_BORMAIO
#include <meteoio/plugins/BormaIO.h>
#endif

#ifdef PLUGIN_CNRMIO
#include <meteoio/plugins/CNRMIO.h>
#endif

#ifdef PLUGIN_COSMOXMLIO
#include <meteoio/plugins/CosmoXMLIO.h>
#endif

#ifdef PLUGIN_DBO
#include <meteoio/plugins/DBO.h>
#endif

#ifdef PLUGIN_IMISIO
#include <meteoio/plugins/ImisIO.h>
#endif

#ifdef PLUGIN_OSHDIO
#include <meteoio/plugins/OshdIO.h>
#endif

#ifdef PLUGIN_GRIBIO
#include <meteoio/plugins/GRIBIO.h>
#endif

#ifdef PLUGIN_GSNIO
#include <meteoio/plugins/GSNIO.h>
#endif

#ifdef PLUGIN_NETCDFIO
#include <meteoio/plugins/NetCDFIO.h>
#endif

#ifdef PLUGIN_PNGIO
#include <meteoio/plugins/PNGIO.h>
#endif

#ifdef PLUGIN_PSQLIO
#include <meteoio/plugins/PSQLIO.h>
#endif

#ifdef PLUGIN_SASEIO
#include <meteoio/plugins/SASEIO.h>
#endif

using namespace std;

namespace mio {
 /**
 * @page data_sources Data input overview
 * The data access is handled by a system of plugins. They all offer the same interface, meaning that a plugin can transparently be replaced by another one. Since they 
 * might rely on third party libraries for accessing the data, they have been created as plugins, that is they are only compiled if requested when configuring the
 * compilation with cmake. A plugin can therefore fail to run if it has not been compiled.
 *
 * Please have a look at the support for \subpage coords "coordinate systems".
 * 
 * @section available_categories Data sources categories
 * Several data sources categories have been defined that can be provided by a different plugin. Each data source category is defined by a specific key in the configuration file (usually, io.ini):
 * - METEO, for meteorological time series
 * - DEM, for Digital Elevation Maps
 * - LANDUSE, for land cover information
 * - GRID2D, for generic 2D grids (they can contain meteo fields and be recognized as such or arbitrary gridded data)
 * - POI, for a list of Points Of Interest that can be used for providing extra information at some specific location (extracting time series at a few selected points, etc)
 *
 * A plugin is "connected" to a given data source category simply by giving its keyword as value for the data source key:
 * @code
 * METEO = SMET
 * DEM = ARC
 * @endcode
 * Each plugin might have its own specific options, meaning that it might require its own keywords. Please check in each plugin documentation the supported options and keys (see links below).
 * Moreover, a given plugin might only support a given category for read or write (for example, PNG: there is no easy and safe way to interpret a given color as a given numeric value without knowing its color scale, so reading a png has been disabled).
 * Finally, the plugins usually don't implement all these categories (for example, ArcGIS file format only describes 2D grids, so the ARC plugin will only deal with 2D grids), so please check what a given plugin implements before connecting it to a specific data source category.
 *
 * @subsection available_plugins Available plugins
 * So far the following plugins have been implemented (by keyword for the io.ini key/value config file). Please read the documentation for each plugin in order to know the plugin-specific keywords:
 * <center><table border="1">
 * <tr><th>Plugin keyword</th><th>Provides</th><th>Description</th><th>Extra requirements</th></tr>
 * <tr><td>\subpage alpug "ALPUG"</td><td>meteo</td><td>data files generated by the %ALPUG meteo stations</td><td></td></tr>
 * <tr><td>\subpage a3d "A3D"</td><td>meteo, poi</td><td>original Alpine3D meteo files</td><td></td></tr>
 * <tr><td>\subpage arc "ARC"</td><td>dem, landuse, grid2d</td><td>ESRI/ARC ascii grid files</td><td></td></tr>
 * <tr><td>\subpage arps "ARPS"</td><td>dem, grid2d, grid3d</td><td>ARPS ascii formatted grids</td><td></td></tr>
 * <tr><td>\subpage borma "BORMA"</td><td>meteo</td><td>Borma xml meteo files</td><td><A HREF="http://libxmlplusplus.sourceforge.net/">libxml++</A></td></tr>
 * <tr><td>\subpage cnrm "CNRM"</td><td>dem, grid2d, meteo</td><td>NetCDF meteorological timeseries following the <A HREF="http://www.cnrm.meteo.fr/?lang=en">CNRM</A> schema</td><td><A HREF="http://www.unidata.ucar.edu/downloads/netcdf/index.jsp">NetCDF-C library</A></td></tr>
 * <tr><td>\subpage cosmoxml "COSMOXML"</td><td>meteo</td><td>MeteoSwiss COSMO's postprocessing XML format</td><td><A HREF="http://xmlsoft.org/">libxml2</A></td></tr>
 * <tr><td>\subpage dbo "DBO"</td><td>meteo</td><td>connects to SLF's DBO web service interface</td><td><A HREF="http://curl.haxx.se/libcurl/">libcurl</A></td></tr>
 * <tr><td>\subpage geotop "GEOTOP"</td><td>meteo</td><td>GeoTop meteo files</td><td></td></tr>
 * <tr><td>\subpage grass "GRASS"</td><td>dem, landuse, grid2d</td><td>Grass grid files</td><td></td></tr>
 * <tr><td>\subpage gribio "GRIB"</td><td>meteo, dem, grid2d</td><td>GRIB meteo grid files</td><td><A HREF="http://www.ecmwf.int/products/data/software/grib_api.html">grib-api</A></td></tr>
 * <tr><td>\subpage gsn "GSN"</td><td>meteo</td><td>connects to the Global Sensor Network web service interface</td><td><A HREF="http://curl.haxx.se/libcurl/">libcurl</A></td></tr>
 * <tr><td>\subpage imis "IMIS"</td><td>meteo</td><td>connects to the IMIS database</td><td><A HREF="http://docs.oracle.com/cd/B12037_01/appdev.101/b10778/introduction.htm">Oracle's OCCI library</A></td></tr>
 * <tr><td>\subpage netcdf "NETCDF"</td><td>dem, grid2d</td><td>NetCDF grids</td><td><A HREF="http://www.unidata.ucar.edu/downloads/netcdf/index.jsp">NetCDF-C library</A></td></tr>
 * <tr><td>\subpage oshd "OSHD"</td><td>meteo</td><td>OSHD generated binary Matlab files</td><td><A HREF="https://sourceforge.net/projects/matio">libmatio</A></td></tr>
 * <tr><td>\subpage pgmio "PGM"</td><td>dem, grid2d</td><td>PGM grid files</td><td></td></tr>
 * <tr><td>\subpage pngio "PNG"</td><td>dem, grid2d</td><td>PNG grid files</td><td><A HREF="http://www.libpng.org/pub/png/libpng.html">libpng</A></td></tr>
 * <tr><td>\subpage psqlio "PSQL"</td><td>meteo</td><td>connects to PostgreSQL database</td><td><A HREF="http://www.postgresql.org/">PostgreSQL</A>'s libpq</td></tr>
 * <tr><td>\subpage sase "SASE"</td><td>meteo</td><td>connects to the SASE database</td><td><A HREF="https://dev.mysql.com/doc/refman/5.0/en/c-api.html">MySQL's C API</A></td></tr>
 * <tr><td>\subpage smetio "SMET"</td><td>meteo, poi</td><td>SMET data files</td><td></td></tr>
 * <tr><td>\subpage snowpack "SNOWPACK"</td><td>meteo</td><td>original SNOWPACK meteo files</td><td></td></tr>
 * </table></center>
 *
 * In order to optimize the data retrieval, the raw data is buffered. This means that up to BUFFER_SIZE days of data will be read at once by the plugin
 * so subsequent reads will not have to get back to the data source (this key is in the [General] section). It is usually a good idea to configure BUFFER_SIZE
 * to the intended duration of the simulation (in days).
 *
 * @section data_manipulations Raw data editing
 * Before any filters, resampling algorithms or data generators are applied, it is possible to edit the original data:
 *     -# \ref data_move "rename certain parameters for all stations;"
 *     -# \ref data_exclusion "exclude/keep certain parameters on a per station basis;"
 *     -# \ref data_merging "merge stations together;"
 *     -# \ref data_copy "make a copy of a certain parameter under a new parameter name for all stations;"
 *     -# \ref data_creation "create certain parameters based on some parametrizations."
 * 
 * @note Please note that the processing order is the following: the MOVE directives are processed first, then the EXCLUDE directives, 
 * then the KEEP directives, then the MERGE directives and finally the COPY directives. The CREATE directives only come after all the raw data
 * has been edited.
 *
 * @subsection data_move 1. Data renaming (MOVE)
 * It is possible to rename a meteorological parameter thanks to the MOVE key. This key can take multiple source names that will be processed in the
 * order of declaration. The syntax is new_name::MOVE = {*space delimited list of original names*}. Original names that are not found in the current
 * dataset will silently be ignored, so it is safe to provide a list that contain many possible names:
 * @code
 * TA::MOVE = air_temp air_temperature temperature_air
 * @endcode
 * This can be used to rename non-standard parameter names into standard ones.
 * 
 * @subsection data_exclusion 2. Data exclusion (EXCLUDE/KEEP)
 * It is possible to exclude specific parameters from given stations (on a per station basis). This is either done by using the station ID (or the '*' wildcard) 
 * followed by "::exclude" as key with a space delimited list of \ref meteoparam "meteorological parameters" to exclude for the station as key.
 * Another possibility is to provide a file containing one station ID per line followed by a space delimited list of \ref meteoparam "meteorological parameters"
 * to exclude for the station (the path to the file can be a relative path and will be properly resolved).
 * 
 * The exact opposite can also be done, excluding ALL parameters except the ones declared with the "::keep" statement (or a file containing one station ID
 * per line followed by a space delimited list of \ref meteoparam "meteorological parameters" to keep for the station).
 *
 * @code
 * WFJ2::EXCLUDE = HS PSUM                       ;inline declaration of parameters exclusion
 * KLO3::KEEP = TA RH VW DW                      ;inline declaration of parameters to keep
 *
 * EXCLUDE_FILE = ../input/meteo/excludes.csv    ;parameters exclusions defined in a separate file
 * KEEP_FILE = ../input/meteo/keeps.csv          ;parameters to keep defined in a separate file
 * @endcode
 *
 * In the second example (relying on a separate file), the file "../input/meteo/excludes.csv" could look like this:
 * @code
 * WFJ2 TA RH
 * KLO3 HS PSUM
 * @endcode
 * 
 * Another example relying on wildcards (the kept/excluded parameters lists are additive):
 * @code
 * *::KEEP = TA RH                               ;all stations will keep TA and RH and reject the other parameters
 * WFJ2::KEEP = HS PSUM                          ;WFJ2 will keep TA and RH as defined above but also HS and PSUM
 * @endcode
 *
 * @subsection data_merging 3. Data merging (MERGE)
 * It is possible to merge different data sets together, with a syntax similar to the Exclude/Keep syntax. This merging occurs <b>after</b> any 
 * EXCLUDE/KEEP commands. This is useful, for example, to provide measurements from different stations that actually share the 
 * same measurement location or to build "composite" station from multiple real stations (in this case, using EXCLUDE and/or KEEP 
 * commands to fine tune how the composite station(s) is/are built). 
 * Please note that the order of declaration defines the priority (ie the first station that has a value for a given parameter has priority). Please also
 * note that only common timestamps will be merged! (ie if the stations have different sampling rates, it might end up that no merge gets performed)
 * 
 * @code
 * STATION1 = STB
 * STATION2 = WFJ2
 * STATION3 = WFJ1
 * STATION4 = DAV1
 * [...]
 * 
 * STB::EXCLUDE = ILWR PSUM
 * WFJ2::KEEP = PSUM ILWR RSWR
 * 
 * STB::MERGE = WFJ2 WFJ1
 * DAV1::MERGE = WFJ2
 * @endcode
 * In order to avoid circular dependencies, a station can NOT receive data from a station AND contribute data to another station. Otherwise, a 
 * station can be merged into multiple other stations. Moreover, the merging strategy can be controlled by setting the MERGE_STRATEGY key in
 * the [Input] section (by default it is "STRICT_MERGE", see MeteoData::Merge_Type).
 * 
 * @note One limitation when handling "extra" parameters (ie parameters that are not in the default \ref meteoparam) is that these extra 
 * parameters must be known from the begining. So if station2 appears later in time with extra parameters, make sure that the buffer size 
 * is large enough to reach all the way to this new station (by setting General::BUFFER_SIZE at least to the number of days from
 * the start of the first station to the start of the second station)
 * 
 * @subsection data_copy 4. Data copy (COPY)
 * It is also possible to duplicate a meteorological parameter as another meteorological parameter. This is done by specifying a COPY key, following the syntax
 * new_name::COPY = existing_parameter. For example:
 * @code
 * VW_avg::COPY = VW
 * @endcode
 * This creates a new parameter VW_avg that starts as an exact copy of the raw data of VW, for each station. This newly created parameter is
 * then processed as any other meteorological parameter (thus going through filtering, generic processing, spatial interpolations). This only current
 * limitation is that the parameter providing the raw data must be defined for all stations (even if filled with nodata, this is good enough).
 *
 * @subsection data_creation 5. Data creation (CREATE)
 * Finally, it is possible to create new data based on some parametrizations. If the requested parameter does not exists, it will be created. Otherwise,
 * any pre-existing data is kept and only missing values in the original data set are filled with the generated values, keeping the original sampling rate. As
 * with all raw data editing, this takes place *before* any filtering/resampling/data generators. As the available algorithms are the same as for the
 * data generators, they are listed in the \ref generators_keywords "data generators section" (but the data creators must be declared in the [Input] section).
 * @code
 * [Input]
 * P::create = STD_PRESS			#the pressure is filled with STD_PRESS if no measured values are available
 * ISWR_POT::create = clearSky_SW		#a new parameter "ISWR_POT" is created and filled with Clear Sky values
 * @endcode
 *
 * @section virtual_stations_section Virtual stations
 * It is possible to use spatially interpolated meteorological fields or time series of 2D grids to extract meteorological time series for a set of points.
 * This is handled as "virtual stations" since the data will seem to originate from points where no station is present. This is described in the 
 * \subpage virtual_stations "virtual stations" page.
 * 
 */

IOInterface* IOHandler::getPlugin(const std::string& plugin_name) const
{
#ifdef PLUGIN_ALPUG
	if (plugin_name == "ALPUG") return new ALPUG(cfg);
#endif
#ifdef PLUGIN_ARCIO
	if (plugin_name == "ARC") return new ARCIO(cfg);
#endif
#ifdef PLUGIN_A3DIO
	if (plugin_name == "A3D") return new A3DIO(cfg);
#endif
#ifdef PLUGIN_ARPSIO
	if (plugin_name == "ARPS") return new ARPSIO(cfg);
#endif
#ifdef PLUGIN_GRASSIO
	if (plugin_name == "GRASS") return new GrassIO(cfg);
#endif
#ifdef PLUGIN_GEOTOPIO
	if (plugin_name == "GEOTOP") return new GeotopIO(cfg);
#endif
#ifdef PLUGIN_SMETIO
	if (plugin_name == "SMET") return new SMETIO(cfg);
#endif
#ifdef PLUGIN_SNIO
	if (plugin_name == "SNOWPACK") return new SNIO(cfg);
#endif
#ifdef PLUGIN_PGMIO
	if (plugin_name == "PGM") return new PGMIO(cfg);
#endif
#ifdef PLUGIN_IMISIO
	if (plugin_name == "IMIS") return new ImisIO(cfg);
#endif
#ifdef PLUGIN_OSHDIO
	if (plugin_name == "OSHD") return new OshdIO(cfg);
#endif
#ifdef PLUGIN_GRIBIO
	if (plugin_name == "GRIB") return new GRIBIO(cfg);
#endif
#ifdef PLUGIN_PNGIO
	if (plugin_name == "PNG") return new PNGIO(cfg);
#endif
#ifdef PLUGIN_BORMAIO
	if (plugin_name == "BORMA") return new BormaIO(cfg);
#endif
#ifdef PLUGIN_COSMOXMLIO
	if (plugin_name == "COSMOXML") return new CosmoXMLIO(cfg);
#endif
#ifdef PLUGIN_DBO
	if (plugin_name == "DBO") return new DBO(cfg);
#endif
#ifdef PLUGIN_CNRMIO
	if (plugin_name == "CNRM") return new CNRMIO(cfg);
#endif
#ifdef PLUGIN_GSNIO
	if (plugin_name == "GSN") return new GSNIO(cfg);
#endif
#ifdef PLUGIN_NETCDFIO
	if (plugin_name == "NETCDF") return new NetCDFIO(cfg);
#endif
#ifdef PLUGIN_PSQLIO
	if (plugin_name == "PSQL") return new PSQLIO(cfg);
#endif
#ifdef PLUGIN_SASEIO
	if (plugin_name == "SASE") return new SASEIO(cfg);
#endif

	return NULL; //no plugin found
}

//this is actually an object factory
IOInterface* IOHandler::getPlugin(const std::string& cfgkey, const std::string& cfgsection)
{
	const std::string op_src = cfg.get(cfgkey, cfgsection);

	if (mapPlugins.find(op_src) == mapPlugins.end()) {
		IOInterface *ioPtr = getPlugin(op_src);
		if (ioPtr==NULL)
			throw IOException("Cannot find plugin " + op_src + " as requested in file " + cfg.getSourceName() + ". Has it been activated through ccmake? Is it declared in IOHandler::getPlugin?", AT);
		else
			mapPlugins[op_src] = ioPtr;
	}

	return mapPlugins[op_src];
}

//Copy constructor
IOHandler::IOHandler(const IOHandler& aio)
           : IOInterface(), cfg(aio.cfg), dataCreator(aio.cfg), mapPlugins(aio.mapPlugins), excluded_params(aio.excluded_params), kept_params(aio.kept_params),
             merge_commands(aio.merge_commands), copy_commands(aio.copy_commands), move_commands(aio.move_commands),
             merged_stations(aio.merged_stations), merge_strategy(aio.merge_strategy), 
             copy_ready(aio.copy_ready), move_ready(aio.move_ready), excludes_ready(aio.excludes_ready), keeps_ready(aio.keeps_ready), merge_ready(aio.merge_ready)
{}

IOHandler::IOHandler(const Config& cfgreader)
           : IOInterface(), cfg(cfgreader), dataCreator(cfgreader), mapPlugins(), excluded_params(), kept_params(),
             merge_commands(), copy_commands(), move_commands(), 
             merged_stations(), merge_strategy(MeteoData::STRICT_MERGE), 
             copy_ready(false), move_ready(false), excludes_ready(false), keeps_ready(false), merge_ready(false)
{
	const std::string merge_strategy_str = cfg.get("MERGE_STRATEGY", "Input", IOUtils::nothrow);
	if (!merge_strategy_str.empty())
		merge_strategy = MeteoData::getMergeType(merge_strategy_str);
}

IOHandler::~IOHandler() throw()
{
	// Get rid of the objects
	std::map<std::string, IOInterface*>::iterator mapit;
	for (mapit = mapPlugins.begin(); mapit!=mapPlugins.end(); ++mapit) {
		delete mapit->second;
	}
}

IOHandler& IOHandler::operator=(const IOHandler& source) {
	if (this != &source) {
		dataCreator = source.dataCreator;
		mapPlugins = source.mapPlugins;
		excluded_params = source.excluded_params;
		kept_params = source.kept_params;
		merge_commands = source.merge_commands;
		merged_stations = source.merged_stations;
		copy_commands = source.copy_commands;
		move_commands = source.move_commands;
		merge_strategy = source.merge_strategy;
		copy_ready = source.copy_ready;
		move_ready = source.move_ready;
		excludes_ready = source.excludes_ready;
		keeps_ready = source.keeps_ready;
		merge_ready = source.merge_ready;
	}
	return *this;
}

void IOHandler::read2DGrid(Grid2DObject& grid_out, const std::string& i_filename)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	plugin->read2DGrid(grid_out, i_filename);
}

void IOHandler::read2DGrid(Grid2DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID2D", "Input");
	plugin->read2DGrid(grid_out, parameter, date);
}

void IOHandler::read3DGrid(Grid3DObject& grid_out, const std::string& i_filename)
{
	IOInterface *plugin = getPlugin("GRID3D", "Input");
	plugin->read3DGrid(grid_out, i_filename);
}

void IOHandler::read3DGrid(Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID3D", "Input");
	plugin->read3DGrid(grid_out, parameter, date);
}

void IOHandler::readDEM(DEMObject& dem_out)
{
	IOInterface *plugin = getPlugin("DEM", "Input");
	plugin->readDEM(dem_out);
	dem_out.update();
}

void IOHandler::readLanduse(Grid2DObject& landuse_out)
{
	IOInterface *plugin = getPlugin("LANDUSE", "Input");
	plugin->readLanduse(landuse_out);
}

void IOHandler::readStationData(const Date& date, STATIONS_SET& vecStation)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readStationData(date, vecStation);
	
	if (!merge_ready) create_merge_map(); 
	merge_stations(vecStation);
}

void IOHandler::readMeteoData(const Date& dateStart, const Date& dateEnd,
                              std::vector<METEO_SET>& vecMeteo)
{
	IOInterface *plugin = getPlugin("METEO", "Input");
	plugin->readMeteoData(dateStart, dateEnd, vecMeteo);

	checkTimestamps(vecMeteo);
	
	if (!move_ready) create_move_map();
	move_params(vecMeteo);
	
	if (!excludes_ready) create_exclude_map();
	exclude_params(vecMeteo);
	
	if (!keeps_ready) create_keep_map();
	keep_params(vecMeteo);
	
	if (!merge_ready) create_merge_map(); 
	merge_stations(vecMeteo);
	
	if (!copy_ready) create_copy_map();
	copy_params(vecMeteo);

	dataCreator.createParameters(vecMeteo);
}

void IOHandler::writeMeteoData(const std::vector<METEO_SET>& vecMeteo,
                               const std::string& name)
{
	IOInterface *plugin = getPlugin("METEO", "Output");
	plugin->writeMeteoData(vecMeteo, name);
}

void IOHandler::readAssimilationData(const Date& date_in, Grid2DObject& da_out)
{
	IOInterface *plugin = getPlugin("DA", "Input");
	plugin->readAssimilationData(date_in, da_out);
}

void IOHandler::readPOI(std::vector<Coords>& pts) {
	IOInterface *plugin = getPlugin("POI", "Input");
	plugin->readPOI(pts);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const std::string& name)
{
	IOInterface *plugin = getPlugin("GRID2D", "Output");
	plugin->write2DGrid(grid_in, name);
}

void IOHandler::write2DGrid(const Grid2DObject& grid_in, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID2D", "Output");
	plugin->write2DGrid(grid_in, parameter, date);
}

void IOHandler::write3DGrid(const Grid3DObject& grid_out, const std::string& options)
{
	IOInterface *plugin = getPlugin("GRID3D", "Output");
	plugin->write3DGrid(grid_out, options);
}

void IOHandler::write3DGrid(const Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	IOInterface *plugin = getPlugin("GRID3D", "Output");
	plugin->write3DGrid(grid_out, parameter, date);
}

/** 
 * @brief check that timestamps are unique and in increasing order
 * @param[in] vecVecMeteo all the data for all the stations
*/
void IOHandler::checkTimestamps(const std::vector<METEO_SET>& vecVecMeteo)
{
	for (size_t stat_idx=0; stat_idx<vecVecMeteo.size(); ++stat_idx) { //for each station
		const size_t nr_timestamps = vecVecMeteo[stat_idx].size();
		if (nr_timestamps==0) continue;

		Date previous_date( vecVecMeteo[stat_idx].front().date );
		for (size_t ii=1; ii<nr_timestamps; ++ii) {
			const Date current_date( vecVecMeteo[stat_idx][ii].date );
			if (current_date<=previous_date) {
				const StationData& station( vecVecMeteo[stat_idx][ii].meta );
				if (current_date==previous_date)
					throw IOException("Error for station \""+station.stationName+"\" ("+station.stationID+") at time "+current_date.toString(Date::ISO)+": timestamps must be unique!", AT);
				else
					throw IOException("Error for station \""+station.stationName+"\" ("+station.stationID+"): jumping from "+previous_date.toString(Date::ISO)+" to "+current_date.toString(Date::ISO), AT);
			}
			previous_date = current_date;
		}
	}
}

void IOHandler::create_merge_map()
{
	merge_ready = true;
	
	std::vector<std::string> merge_keys;
	const size_t nrOfStations = cfg.findKeys(merge_keys, "::MERGE", "Input", true);
	for (size_t ii=0; ii<nrOfStations; ++ii) {
		const size_t found = merge_keys[ii].find_first_of(":");
		if (found==std::string::npos) continue;

		const std::string station( IOUtils::strToUpper(merge_keys[ii].substr(0,found)) );
		std::vector<std::string> vecString;
		cfg.getValue(merge_keys[ii], "Input", vecString);
		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+merge_keys[ii]+"\"", AT);
		
		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
			const std::vector<std::string>::const_iterator vec_it = find (merged_stations.begin(), merged_stations.end(), *it);
			if (vec_it==merged_stations.end()) merged_stations.push_back( *it ); //this station will be merged into another one
		}
		merge_commands[ station ] = vecString;
	}
	
	//sort the merged_stations vector so searches will be faster
	std::sort(merged_stations.begin(), merged_stations.end());
	
	//make sure there is no "chain merge": station A merging station B and station C merging station A
	std::map< std::string, std::vector<std::string> >::iterator it_dest;
	for(it_dest=merge_commands.begin(); it_dest!=merge_commands.end(); ++it_dest) {
		const std::string stationID( it_dest->first );
		if (std::binary_search(merged_stations.begin(), merged_stations.end(), stationID))
			throw InvalidArgumentException("\'chain merge\' detected for station \'"+stationID+"\', this is not supported (see documentation)", AT);
	}
}

//merge stations that have identical names
void IOHandler::merge_stations(STATIONS_SET& vecStation) const
{
	if (merge_commands.empty()) return;
	
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		const std::string toStationID( IOUtils::strToUpper( vecStation[ii].stationID ) );
		//we do not support "chain merge": station A merging station B and station C merging station A
		if ( std::find(merged_stations.begin(), merged_stations.end(), toStationID)!=merged_stations.end() ) continue;
		
		const std::map< string, vector<string> >::const_iterator it = merge_commands.find( toStationID );
		if (it == merge_commands.end()) continue; //no merge commands for this station

		const std::vector<std::string> merge_from( it->second );
		for (std::vector<std::string>::const_iterator it_set=merge_from.begin(); it_set != merge_from.end(); ++it_set) {
			const std::string fromStationID( *it_set );
			
			bool found = false;
			for (size_t jj=0; jj<vecStation.size(); jj++) {
				const std::string curr_station( IOUtils::strToUpper(vecStation[jj].stationID) );
				if (curr_station==fromStationID) {
					vecStation[ii].merge( vecStation[jj] );
					found = true;
				}
			}
			if (!found)
				throw InvalidArgumentException("Station ID '"+fromStationID+"' not found when merging toward station '"+toStationID+"'. Consider increasing BUFFER_SIZE!", AT);
		}
	}
	
	//remove the stations that have been merged into other ones
	for (size_t ii=0; ii<vecStation.size(); ii++) {
		const std::string stationID( IOUtils::strToUpper( vecStation[ii].stationID ) );
		const std::vector<std::string>::const_iterator it = std::find(merged_stations.begin(), merged_stations.end(), stationID);
		if ( it!=merged_stations.end() ) {
			std::swap( vecStation[ii], vecStation.back() );
			vecStation.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
}

//in this implementation, we consider that the station name does NOT change over time
void IOHandler::merge_stations(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (merge_commands.empty()) return;
	
	for (size_t ii=0; ii<vecVecMeteo.size(); ii++) { //loop over the stations
		if (vecVecMeteo[ii].empty())  continue;
		const std::string toStationID( IOUtils::strToUpper(vecVecMeteo[ii][0].meta.stationID) );
		//we do not support "chain merge": station A merging station B and station C merging station A
		if ( std::find(merged_stations.begin(), merged_stations.end(), toStationID)!=merged_stations.end() ) continue;
		
		const std::map< std::string, std::vector<std::string> >::const_iterator it = merge_commands.find( toStationID );
		if (it == merge_commands.end()) continue; //no merge commands for this station

		const std::vector<std::string> merge_from( it->second );
		for (std::vector<std::string>::const_iterator it_set=merge_from.begin(); it_set != merge_from.end(); ++it_set) {
			const std::string fromStationID( *it_set );
			
			bool found = false;
			for (size_t jj=0; jj<vecVecMeteo.size(); jj++) { //loop over the available stations in the current dataset
				if (vecVecMeteo[jj].empty()) continue;
				const std::string curr_station( IOUtils::strToUpper(vecVecMeteo[jj][0].meta.stationID) );
				if (curr_station==fromStationID) {
					MeteoData::mergeTimeSeries(vecVecMeteo[ii], vecVecMeteo[jj], static_cast<MeteoData::Merge_Type>(merge_strategy)); //merge timeseries for the two stations
					found = true;
				}
			}
			if (!found)
				throw InvalidArgumentException("Station ID '"+fromStationID+"' not found when merging toward station '"+toStationID+"'. Consider increasing BUFFER_SIZE!", AT);
		}
	}
	
	//remove the stations that have been merged into other ones
	for (size_t ii=0; ii<vecVecMeteo.size(); ii++) {
		if (vecVecMeteo[ii].empty())  continue;
		const std::string stationID( IOUtils::strToUpper(vecVecMeteo[ii][0].meta.stationID) );
		const std::vector<std::string>::const_iterator it = std::find(merged_stations.begin(), merged_stations.end(), stationID);
		if ( it!=merged_stations.end() ) {
			std::swap( vecVecMeteo[ii], vecVecMeteo.back() );
			vecVecMeteo.pop_back();
			ii--; //in case we have multiple identical stations ID
		}
	}
}

void IOHandler::create_exclude_map()
{
	excludes_ready = true;
	const std::string exclude_file = cfg.get("EXCLUDE_FILE", "Input", IOUtils::nothrow);

	if (!exclude_file.empty()) {
		//if this is a relative path, prefix the path with the current path
		const std::string prefix = ( FileUtils::isAbsolutePath(exclude_file) )? "" : cfg.getConfigRootDir()+"/";
		const std::string path( FileUtils::getPath(prefix+exclude_file, true) );  //clean & resolve path
		const std::string filename( path + "/" + FileUtils::getFilename(exclude_file) );

		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		std::ifstream fin(filename.c_str(), std::ifstream::in);
		if (fin.fail()) throw AccessException(filename, AT);

		try {
			const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

			std::vector<std::string> tmpvec;
			std::string line;

			while (!fin.eof()) { //Go through file
				getline(fin, line, eoln); //read complete line meta information
				IOUtils::stripComments(line);
				const size_t ncols = IOUtils::readLineToVec(line, tmpvec, ' ');

				if (ncols > 1) {
					for (std::vector<std::string>::iterator it = tmpvec.begin()+1; it != tmpvec.end(); ++it) {
						IOUtils::toUpper( *it );
					}

					const std::set<std::string> tmpset(tmpvec.begin()+1, tmpvec.end());
					excluded_params[ IOUtils::strToUpper(tmpvec[0]) ] = tmpset;
				}
			}
		} catch (const std::exception&) {
			fin.close();
			throw;
		}

		fin.close();
	}

	std::vector<std::string> exclude_keys;
	const size_t nrOfStations = cfg.findKeys(exclude_keys, "::EXCLUDE", "Input", true);
	for (size_t ii=0; ii<nrOfStations; ++ii) {
		const size_t found = exclude_keys[ii].find_first_of(":");
		if (found==std::string::npos) continue;

		const std::string station( IOUtils::strToUpper(exclude_keys[ii].substr(0,found)) );
		std::vector<std::string> vecString;
		cfg.getValue(exclude_keys[ii], "Input", vecString);
		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+exclude_keys[ii]+"\"", AT);
		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
		}

		const std::set<std::string> tmpset(vecString.begin(), vecString.end());
		excluded_params[ station ] = tmpset;
	}
	
	//Handle "*" wildcard: add the params to all other declared stations
	std::map< std::string, std::set<std::string> >::const_iterator it_station = excluded_params.find("*");
	if (it_station!=excluded_params.end()) {
		const std::set<std::string> wildcard( excluded_params["*"] );
		for (it_station=excluded_params.begin(); it_station!=excluded_params.end(); ++it_station) {
			std::set<std::string> params( it_station->second );
			
			for (std::set<std::string>::iterator it=wildcard.begin(); it!=wildcard.end(); ++it)
				params.insert( *it ); //merging: keep in mind that a set can not contain duplicates
			
			excluded_params[ it_station->first ] = params;
		}
	}
}

void IOHandler::create_keep_map()
{
	keeps_ready = true;
	const std::string keep_file = cfg.get("KEEP_FILE", "Input", IOUtils::nothrow);

	if (!keep_file.empty()) {
		//if this is a relative path, prefix the path with the current path
		const std::string prefix = ( FileUtils::isAbsolutePath(keep_file) )? "" : cfg.getConfigRootDir()+"/";
		const std::string path( FileUtils::getPath(prefix+keep_file, true) );  //clean & resolve path
		const std::string filename( path + "/" + FileUtils::getFilename(keep_file) );

		if (!FileUtils::fileExists(filename)) throw AccessException(filename, AT); //prevent invalid filenames
		std::ifstream fin(filename.c_str(), std::ifstream::in);
		if (fin.fail()) throw AccessException(filename, AT);

		try {
			const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

			std::vector<std::string> tmpvec;
			std::string line;

			while (!fin.eof()) { //Go through file
				getline(fin, line, eoln); //read complete line meta information
				IOUtils::stripComments(line);
				const size_t ncols = IOUtils::readLineToVec(line, tmpvec, ' ');

				if (ncols > 1) {
					for (vector<string>::iterator it = tmpvec.begin()+1; it != tmpvec.end(); ++it) {
						IOUtils::toUpper( *it );
					}

					const set<string> tmpset(tmpvec.begin()+1, tmpvec.end());
					kept_params[ IOUtils::strToUpper(tmpvec[0]) ] = tmpset;
				}
			}
		} catch (const std::exception&) {
			fin.close();
			throw;
		}

		fin.close();
	}

	std::vector<std::string> keep_keys;
	const size_t nrOfStations = cfg.findKeys(keep_keys, "::KEEP", "Input", true);
	for (size_t ii=0; ii<nrOfStations; ++ii) {
		const size_t found = keep_keys[ii].find_first_of(":");
		if (found==std::string::npos) continue;

		const std::string station( IOUtils::strToUpper(keep_keys[ii].substr(0,found)) );
		std::vector<std::string> vecString;
		cfg.getValue(keep_keys[ii], "Input", vecString);
		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+keep_keys[ii]+"\"", AT);
		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
		}

		const std::set<std::string> tmpset(vecString.begin(), vecString.end());
		kept_params[ station ] = tmpset;
	}
	
	//Handle "*" wildcard: add the params to all other declared stations
	std::map< std::string, std::set<std::string> >::const_iterator it_station = kept_params.find("*");
	if (it_station!=kept_params.end()) {
		const std::set<std::string> wildcard( kept_params["*"] );
		for (it_station=kept_params.begin(); it_station!=kept_params.end(); ++it_station) {
			std::set<std::string> params( it_station->second );
			
			for (std::set<std::string>::iterator it=wildcard.begin(); it!=wildcard.end(); ++it)
				params.insert( *it ); //merging: keep in mind that a set can not contain duplicates
			
			kept_params[ it_station->first ] = params;
		}
	}
}

/**
* @brief reset to nodata the parameters marked as EXCLUDE on a per station basis
*/
void IOHandler::exclude_params(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (excluded_params.empty()) return;

	for (size_t station=0; station<vecVecMeteo.size(); ++station) { //loop over the stations
		if (vecVecMeteo[station].empty()) continue;
		const std::string stationID( IOUtils::strToUpper(vecVecMeteo[station][0].meta.stationID) );
		std::map< std::string, std::set<std::string> >::const_iterator it = excluded_params.find(stationID);
		if (it == excluded_params.end()) {
			it = excluded_params.find("*"); //fallback: is there a wildcard like "*::KEEP"?
			if (it == excluded_params.end()) continue;
		}

		const std::set<std::string> excluded( it->second );

		for (size_t ii=0; ii<vecVecMeteo[station].size(); ++ii) { //loop over the timesteps
			for (std::set<std::string>::const_iterator it_set=excluded.begin(); it_set != excluded.end(); ++it_set) {
				const std::string param( *it_set );
				if (vecVecMeteo[station][ii].param_exists(param))
					vecVecMeteo[station][ii](param) = IOUtils::nodata;
			}
		}
	}
}

/**
* @brief only keep the parameters marked as KEEP on a per station basis
*/
void IOHandler::keep_params(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (kept_params.empty()) return;

	for (size_t station=0; station<vecVecMeteo.size(); ++station) { //loop over the stations
		if (vecVecMeteo[station].empty()) continue;
		
		const std::string stationID( IOUtils::strToUpper(vecVecMeteo[station][0].meta.stationID) );
		std::map< std::string, std::set<std::string> >::const_iterator it = kept_params.find(stationID);
		if (it == kept_params.end()) {
			it = kept_params.find("*"); //fallback: is there a wildcard like "*::KEEP"?
			if (it == kept_params.end()) continue;
		}

		const std::set<std::string> kept( it->second );
		
		for (size_t ii=0; ii<vecVecMeteo[station].size(); ++ii) {
			MeteoData& md_ref( vecVecMeteo[station][ii] );
			MeteoData md( md_ref );
			md.reset(); //delete all meteo fields
			
			for (std::set<std::string>::const_iterator it_set=kept.begin(); it_set != kept.end(); ++it_set) { //loop over the parameters to keep
				const std::string param( *it_set);
				if (!md.param_exists(param)) continue;
				 md(param) = md_ref(param);
			}
			
			//copy back the new object into vecVecMeteo
			md_ref = md;
		}
	}
}

/**
* Parse [Input] section for potential parameters that the user wants
* renamed (as '%%::MOVE = %%')
*/
void IOHandler::create_move_map()
{
	std::vector<std::string> move_keys;
	const size_t nrOfMatches = cfg.findKeys(move_keys, "::MOVE", "Input", true); //search anywhere in key

	for (size_t ii=0; ii<nrOfMatches; ++ii) {
		const std::string dest_param( move_keys[ii].substr( 0, move_keys[ii].find_first_of(":") ) );
		std::vector<std::string> vecString;
		cfg.getValue(move_keys[ii], "Input", vecString); //multiple source can be provided
		
		if (vecString.empty()) throw InvalidArgumentException("Empty value for key \""+move_keys[ii]+"\"", AT);
		for (vector<string>::iterator it = vecString.begin(); it != vecString.end(); ++it) {
			IOUtils::toUpper( *it );
		}
		
		const std::set<std::string> tmpset(vecString.begin(), vecString.end());
		move_commands[ dest_param ] = tmpset;
	}
	
	move_ready = true;
}

/**
* This procedure runs through the MeteoData objects in vecMeteo and according to user
* configuration renames a certain present meteo parameter to another one, named by the
* user in the [Input] section of the io.ini, e.g.
* [Input]
* TA::MOVE = air_temp air_temperature
* means that TA will be the name of a new parameter in MeteoData with the copied value
* of the original parameter air_temp or air_temperature
*/
void IOHandler::move_params(std::vector< METEO_SET >& vecMeteo) const
{
	if (move_commands.empty()) return; //Nothing configured

	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;
		
		std::map< std::string, std::set<std::string> >::const_iterator param = move_commands.begin();
		for (; param!=move_commands.end(); ++param) { //loop over all the MOVE commands
			const std::string dest_param( param->first );
			const std::set<std::string> src( param->second );
			
			for (std::set<std::string>::const_iterator it_set=src.begin(); it_set != src.end(); ++it_set) { //loop over the parameters to move
				const std::string src_param( *it_set );
				const size_t src_index = vecMeteo[station].front().getParameterIndex( src_param );
				if (src_index == IOUtils::npos) continue; //no such parameter for this station, skipping

				for (size_t jj=0; jj<vecMeteo[station].size(); ++jj) {
					const size_t dest_index = vecMeteo[station][jj].addParameter( dest_param ); //either add or just return the proper index
					vecMeteo[station][jj]( dest_index ) = vecMeteo[station][jj]( src_index );
					vecMeteo[station][jj]( src_index ) = IOUtils::nodata; 
				}
			}
		}
	}
}


/**
* Parse [Input] section for potential parameters that the user wants
* duplicated (as '%%::COPY = %%')
*/
void IOHandler::create_copy_map()
{
	std::vector<std::string> copy_keys;
	const size_t nrOfMatches = cfg.findKeys(copy_keys, "::COPY", "Input", true); //search anywhere in key

	for (size_t ii=0; ii<nrOfMatches; ++ii) {
		const std::string dest_param( copy_keys[ii].substr( 0, copy_keys[ii].find_first_of(":") ) );
		const std::string src_param = cfg.get(copy_keys[ii], "Input");
		if (!dest_param.empty() && !src_param.empty())
			copy_commands[ dest_param ] = src_param;
	}
	
	copy_ready = true;
}

/**
* This procedure runs through the MeteoData objects in vecMeteo and according to user
* configuration copies a certain present meteo parameter to another one, named by the
* user in the [Input] section of the io.ini, e.g.
* [Input]
* TA2::COPY = TA
* means that TA2 will be the name of a new parameter in MeteoData with the copied value
* of the meteo parameter MeteoData::TA
*/
void IOHandler::copy_params(std::vector< METEO_SET >& vecMeteo) const
{
	if (copy_commands.empty()) return; //Nothing configured

	for (size_t station=0; station<vecMeteo.size(); ++station) { //for each station
		if (vecMeteo[station].empty()) continue;
		
		std::map< std::string, std::string >::const_iterator param = copy_commands.begin();
		for (; param!=copy_commands.end(); ++param) {
			const size_t src_index = vecMeteo[station].front().getParameterIndex(param->second);
			if (src_index == IOUtils::npos) {
				const std::string stationID( vecMeteo[station].front().meta.stationID );
				throw InvalidArgumentException("Station "+stationID+" has no parameter '"+param->second+"' to copy", AT);
			}
			
			const std::string dest( param->first );
			for (size_t jj=0; jj<vecMeteo[station].size(); ++jj) { //for each MeteoData object of one station
				const size_t dest_index = vecMeteo[station][jj].addParameter( dest );
				vecMeteo[station][jj]( dest_index ) = vecMeteo[station][jj]( src_index );
			}
		}
	}
}

const std::string IOHandler::toString() const
{
	std::ostringstream os;
	os << "<IOHandler>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";

	os << "<mapPlugins>\n";
	std::map<std::string, IOInterface*>::const_iterator it1;
	for (it1=mapPlugins.begin(); it1 != mapPlugins.end(); ++it1) {
		os << setw(10) << it1->first << " = " << hex <<  it1->second << dec << "\n";
	}
	os << "</mapPlugins>\n";

	if (!excluded_params.empty()) {
		os << "<excluded_params>\n";
		std::map< std::string, std::set<std::string> >::const_iterator it_exc;
		for (it_exc=excluded_params.begin(); it_exc != excluded_params.end(); ++it_exc) {
			os << setw(10) << it_exc->first << " = ";
			std::set<std::string>::const_iterator it_set;
			for (it_set=(it_exc->second).begin(); it_set != (it_exc->second).end(); ++it_set)
				os << *it_set << " ";
			os << "\n";
		}
		os << "</excluded_params>\n";
	}
	
	if (!merge_commands.empty()) {
		os << "<merge_commands>\n";
		std::map< std::string, std::vector<std::string> >::const_iterator it_merge;
		for (it_merge=merge_commands.begin(); it_merge != merge_commands.end(); ++it_merge) {
			os << setw(10) << it_merge->first << " <- ";
			for (size_t ii=0; ii<it_merge->second.size(); ++ii)
				os << it_merge->second[ii] << " ";
			os << "\n";
		}
		os << "</merge_commands>\n";
	}

	os << dataCreator.toString();

	os << "</IOHandler>\n";
	return os.str();
}

} //end namespace
