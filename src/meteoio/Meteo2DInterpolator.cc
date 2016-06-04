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

#include <meteoio/Meteo2DInterpolator.h>

using namespace std;

namespace mio {

Meteo2DInterpolator::Meteo2DInterpolator(const Config& i_cfg, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager)
                    : cfg(i_cfg), tsmanager(i_tsmanager), gridsmanager(i_gridsmanager),
                      grid_buffer(0), mapAlgorithms(),
                      v_params(), v_coords(), v_stations(), virtual_point_cache(),
                      algorithms_ready(false), use_full_dem(false), downscaling(false), virtual_stations(false)
{
	size_t max_grids = 10; //default number of grids to keep in buffer
	cfg.getValue("BUFF_GRIDS", "Interpolations2D", max_grids, IOUtils::nothrow);
	grid_buffer.setMaxGrids(max_grids);
	
	setAlgorithms();
	cfg.getValue("Virtual_stations", "Input", virtual_stations, IOUtils::nothrow);
	if (virtual_stations) {
		initVirtualStations();
	}
}

Meteo2DInterpolator::~Meteo2DInterpolator()
{
	std::map<std::string, std::vector<InterpolationAlgorithm*> >::iterator iter;
	for (iter = mapAlgorithms.begin(); iter != mapAlgorithms.end(); ++iter) {
		const vector<InterpolationAlgorithm*>& vecAlgs = iter->second;
		for (size_t ii=0; ii<vecAlgs.size(); ++ii)
			delete vecAlgs[ii];
	}
}

/* By reading the Config object build up a list of user configured algorithms
* for each MeteoData::Parameters parameter (i.e. each member variable of MeteoData like ta, p, psum, ...)
* Concept of this constructor: loop over all MeteoData::Parameters and then look
* for configuration of interpolation algorithms within the Config object.
*/
void Meteo2DInterpolator::setAlgorithms()
{
//HACK set callback to internal iomanager for virtual stations and downsampling!
	set<string> set_of_used_parameters;
	get_parameters(cfg, set_of_used_parameters);

	set<string>::const_iterator it;
	for (it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it) {
		const std::string parname = *it;
		std::vector<std::string> tmpAlgorithms;
		const size_t nrOfAlgorithms = getAlgorithmsForParameter(cfg, parname, tmpAlgorithms);

		std::vector<InterpolationAlgorithm*> vecAlgorithms(nrOfAlgorithms);
		for (size_t jj=0; jj<nrOfAlgorithms; jj++) {
			std::vector<std::string> vecArgs;
			getArgumentsForAlgorithm(parname, tmpAlgorithms[jj], vecArgs);
			vecAlgorithms[jj] = AlgorithmFactory::getAlgorithm( tmpAlgorithms[jj], *this, vecArgs, tsmanager, gridsmanager);
		}

		if (nrOfAlgorithms>0) {
			mapAlgorithms[parname] = vecAlgorithms;
		}
	}
	algorithms_ready = true;
}

//get a list of all meteoparameters referenced in the Interpolations2D section
size_t Meteo2DInterpolator::get_parameters(const Config& cfg, std::set<std::string>& set_parameters)
{
	std::vector<std::string> vec_keys;
	cfg.findKeys(vec_keys, std::string(), "Interpolations2D");

	for (size_t ii=0; ii<vec_keys.size(); ii++) {
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			const string tmp = vec_keys[ii].substr(0,found);
			set_parameters.insert( IOUtils::strToUpper(tmp) );
		}
	}

	return set_parameters.size();
}

void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                                      Grid2DObject& result)
{
	std::string InfoString;
	interpolate(date, dem, meteoparam, result, InfoString);
}

void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                                      Grid2DObject& result, std::string& InfoString)
{
	if (!algorithms_ready)
		setAlgorithms();

	const string param_name = MeteoData::getParameterName(meteoparam);
	
	//Get grid from buffer if it exists
	std::ostringstream grid_hash;
	grid_hash << dem.llcorner.printLatLon() << " " << dem.getNx() << "x" << dem.getNy() << " @" << dem.cellsize << " " << date.toString(Date::ISO) << " " << param_name;
	if (grid_buffer.get(result, grid_hash.str(), InfoString))
		return;

	//Show algorithms to be used for this parameter
	const map<string, vector<InterpolationAlgorithm*> >::iterator it = mapAlgorithms.find(param_name);
	if (it==mapAlgorithms.end()) {
		throw IOException("No interpolation algorithms configured for parameter "+param_name, AT);
	}

	//look for algorithm with the highest quality rating
	const vector<InterpolationAlgorithm*>& vecAlgs = it->second;
	double maxQualityRating = -1.;
	size_t bestalgorithm = 0;
	for (size_t ii=0; ii < vecAlgs.size(); ++ii){
		const double rating = vecAlgs[ii]->getQualityRating(date, meteoparam);
		if ((rating != 0.0) && (rating > maxQualityRating)) {
			//we use ">" so that in case of equality, the first choice will be kept
			bestalgorithm = ii;
			maxQualityRating = rating;
		}
	}

	//finally execute the algorithm with the best quality rating or throw an exception
	if (maxQualityRating<=0.0)
		throw IOException("No interpolation algorithm with quality rating >0 found for parameter "+param_name+" on "+date.toString(Date::ISO_TZ), AT);
	vecAlgs[bestalgorithm]->calculate(dem, result);
	InfoString = vecAlgs[bestalgorithm]->getInfo();

	//Run soft min/max filter for RH, PSUM and HS
	if (meteoparam == MeteoData::RH){
		Meteo2DInterpolator::checkMinMax(0.0, 1.0, result);
	} else if (meteoparam == MeteoData::PSUM){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	} else if (meteoparam == MeteoData::HS){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	} else if (meteoparam == MeteoData::VW){
		Meteo2DInterpolator::checkMinMax(0.0, 10000.0, result);
	}

	//save grid in buffer
	grid_buffer.push(result, grid_hash.str(), InfoString);
}

//HACK make sure that skip_virtual_stations = true before calling this method when using virtual stations!
void Meteo2DInterpolator::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string)
{
	result.clear();
	vector<Coords> vec_coords(in_coords);

	if (use_full_dem) {
		Grid2DObject result_grid;
		interpolate(date, dem, meteoparam, result_grid, info_string);
		const bool gridify_success = dem.gridify(vec_coords);
		if (!gridify_success)
			throw InvalidArgumentException("Coordinate given to interpolate is outside of dem", AT);

		for (size_t ii=0; ii<vec_coords.size(); ii++) {
			//we know the i,j are positive because of gridify_success
			const size_t pt_i = static_cast<size_t>( vec_coords[ii].getGridI() );
			const size_t pt_j = static_cast<size_t>( vec_coords[ii].getGridJ() );
			result.push_back( result_grid(pt_i,pt_j) );
		}
	} else {
		for (size_t ii=0; ii<vec_coords.size(); ii++) {
			const bool gridify_success = dem.gridify(vec_coords[ii]);
			if (!gridify_success)
				throw InvalidArgumentException("Coordinate given to interpolate is outside of dem", AT);

			//we know the i,j are positive because of gridify_success
			const size_t pt_i = static_cast<size_t>( vec_coords[ii].getGridI() );
			const size_t pt_j = static_cast<size_t>( vec_coords[ii].getGridJ() );

			//Make new DEM with just one point, namely the one specified by vec_coord[ii]
			//Copy all other properties of the big DEM into the new one
			DEMObject one_point_dem(dem, pt_i, pt_j, 1, 1, false);

			one_point_dem.min_altitude = dem.min_altitude;
			one_point_dem.max_altitude = dem.max_altitude;
			one_point_dem.min_slope = dem.min_slope;
			one_point_dem.max_slope = dem.max_slope;
			one_point_dem.min_curvature = dem.min_curvature;
			one_point_dem.max_curvature = dem.max_curvature;

			Grid2DObject result_grid;
			interpolate(date, one_point_dem, meteoparam, result_grid, info_string);
			result.push_back(result_grid(0,0));
		}
	}
}

size_t Meteo2DInterpolator::getAlgorithmsForParameter(const Config& cfg, const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	// This function retrieves the user defined interpolation algorithms for
	// parameter 'parname' by querying the Config object
	vecAlgorithms.clear();
	std::vector<std::string> vecKeys;
	cfg.findKeys(vecKeys, parname+"::algorithms", "Interpolations2D");

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::algorithms in config file", AT);;

	if (vecKeys.empty())
		return 0;

	cfg.getValue(vecKeys[0], "Interpolations2D", vecAlgorithms, IOUtils::nothrow);
	return vecAlgorithms.size();
}

size_t Meteo2DInterpolator::getArgumentsForAlgorithm(const std::string& param,
                                                     const std::string& algorithm,
                                                     std::vector<std::string>& vecArgs) const
{
	vecArgs.clear();
	const string keyname = param +"::"+ algorithm;
	cfg.getValue(keyname, "Interpolations2D", vecArgs, IOUtils::nothrow);

	return vecArgs.size();
}

void Meteo2DInterpolator::checkMinMax(const double& minval, const double& maxval, Grid2DObject& gridobj)
{
	const size_t nxy = gridobj.getNx() * gridobj.getNy();

	for (size_t ii=0; ii<nxy; ii++){
		double& value = gridobj(ii);
		if (value == IOUtils::nodata){
			continue;
		}
		if (value < minval) {
			value = minval;
		} else if (value > maxval) {
			value = maxval;
		}
	}
}

void Meteo2DInterpolator::check_projections(const DEMObject& dem, const std::vector<MeteoData>& vec_meteo)
{
	//check that the stations are using the same projection as the dem
	for (size_t ii=0; ii<vec_meteo.size(); ii++) {
		const StationData& meta = vec_meteo[ii].meta;
		if (!meta.position.isSameProj(dem.llcorner)) {
			std::ostringstream os;
			std::string type, args;
			meta.position.getProj(type, args);
			os << "Station " << meta.stationID << " is using projection (" << type << " " << args << ") ";
			dem.llcorner.getProj(type, args);
			os << "while DEM is using projection ("<< type << " " << args << ") ";
			throw IOException(os.str(), AT);
		}
	}
}

//get the stations' data to use for downscaling (=true measurements)
size_t Meteo2DInterpolator::getVirtualMeteoData(const vstations_policy& strategy, const Date& i_date, METEO_SET& vecMeteo)
{
	if (strategy==VSTATIONS) {
		return getVirtualStationsData(i_date, vecMeteo);
	} else if (strategy==DOWNSCALING) {
		//extract all grid points
		return 0; //hack
	} else if (strategy==SMART_DOWNSCALING) {
		//for each parameter:
		//call iomanager.read2DGrid()
		//extract relevant points and fill vecMeteo
		//and loop over all grids!
		//
		//Questions/tricks:
		//   throw exception if grids don't match between parameters
		//   should we recompute which points to take between each time steps? (ie more flexibility)
		//   should the generated virtual stations data be filtered? (useful for applying corrections)
		//          if so: have an own iomanager and use iomanager.push_meteo_data()
		return 0; //hack
	}

	throw UnknownValueException("Unknown virtual station strategy", AT);
}

void Meteo2DInterpolator::initVirtualStations()
{
	if(!cfg.keyExists("DEM", "Input"))
		throw NoAvailableDataException("In order to use virtual stations, please provide a DEM!", AT);
	DEMObject dem;
	gridsmanager.readDEM(dem);

	//get virtual stations coordinates
	std::string coordin, coordinparam, coordout, coordoutparam;
	IOUtils::getProjectionParameters(cfg, coordin, coordinparam, coordout, coordoutparam);

	std::vector<std::string> vecStation;
	cfg.getValues("Vstation", "INPUT", vecStation);
	for(size_t ii=0; ii<vecStation.size(); ii++) {
		//The coordinate specification is given as either: "easting northing epsg" or "lat lon"
		Coords tmp(coordin, coordinparam, vecStation[ii]);
		if(!tmp.isNodata())
			v_coords.push_back( tmp );
	}

	//create stations' metadata
	for(size_t ii=0; ii<v_coords.size(); ii++) {
		if(!dem.gridify(v_coords[ii])) {
			ostringstream ss;
			ss << "Virtual station \"" << vecStation[ii] << "\" is not contained is provided DEM";
			throw NoAvailableDataException(ss.str(), AT);
		}

		const size_t i = v_coords[ii].getGridI(), j = v_coords[ii].getGridJ();
		v_coords[ii].setAltitude(dem(i,j), false);

		ostringstream name;
		name << "Virtual_Station_" << ii+1;
		ostringstream id;
		id << "VIR" << ii+1;
		StationData sd(v_coords[ii], id.str(), name.str());
		sd.setSlope(dem.slope(i,j), dem.azi(i,j));

		v_stations.push_back( sd );
	}

	cfg.getValue("Interpol_Use_Full_DEM", "Input", use_full_dem, IOUtils::nothrow);
}

size_t Meteo2DInterpolator::getVirtualStationsData(const Date& i_date, METEO_SET& vecMeteo)
{ //HACK use own private iomanager and do caching in its private one
	vecMeteo.clear();

	// Check if data is available in cache
	const map<Date, vector<MeteoData> >::const_iterator it = virtual_point_cache.find(i_date);
	if (it != virtual_point_cache.end()){
		vecMeteo = it->second;
		return vecMeteo.size();
	}

	//get data from real input stations
	METEO_SET vecTrueMeteo;
	tsmanager.getMeteoData(i_date, vecTrueMeteo);
	if (vecTrueMeteo.empty()) return 0;

	if (v_params.empty()) {
		//get parameters to interpolate if not already done
		//we need valid data in order to handle extra parameters
		std::vector<std::string> vecStr;
		cfg.getValue("Virtual_parameters", "Input", vecStr);
		for (size_t ii=0; ii<vecStr.size(); ii++) {
			v_params.push_back( vecTrueMeteo[0].getParameterIndex(vecStr[ii]) );
		}
	}

	//create stations without measurements
	for (size_t ii=0; ii<v_stations.size(); ii++) {
		MeteoData md(i_date, v_stations[ii]);
		vecMeteo.push_back( md );
	}

	//fill meteo parameters
	DEMObject dem;
	gridsmanager.readDEM(dem); //this is not a big deal since it will be in the buffer
	string info_string;
	for (size_t param=0; param<v_params.size(); param++) {
		std::vector<double> result;
		interpolate(i_date, dem, static_cast<MeteoData::Parameters>(v_params[param]), v_coords, result, info_string);
		for (size_t ii=0; ii<v_coords.size(); ii++)
			vecMeteo[ii](v_params[param]) = result[ii];
	}

	return vecMeteo.size();
}


const std::string Meteo2DInterpolator::toString() const {
	ostringstream os;
	os << "<Meteo2DInterpolator>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "TimeSeriesManager& tsmanager = "  << hex << &tsmanager << dec << "\n";
	os << "GridsManager& gridsmanager = "  << hex << &gridsmanager << dec << "\n";

	os << "Spatial resampling algorithms:\n";
	std::map<std::string, std::vector<InterpolationAlgorithm*> >::const_iterator iter;
	for (iter = mapAlgorithms.begin(); iter != mapAlgorithms.end(); ++iter) {
		os << setw(10) << iter->first << "::";
		for (size_t jj=0; jj<iter->second.size(); jj++) {
			os << iter->second[jj]->algo << " ";
		}
		os << "\n";
	}

	//cache content
	os << grid_buffer.toString();
	os << "</Meteo2DInterpolator>\n";
	return os.str();
}

} //namespace
