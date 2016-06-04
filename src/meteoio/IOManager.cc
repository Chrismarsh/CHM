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

#include <meteoio/IOManager.h>

using namespace std;

namespace mio {

IOManager::IOManager(const std::string& filename_in) : cfg(filename_in), iohandler(cfg),
                                                       tsmanager(iohandler, cfg), gridsmanager(iohandler, cfg), interpolator(cfg, tsmanager, gridsmanager),
                                                       downscaling(false), virtual_stations(false)
{
	initIOManager();
}

IOManager::IOManager(const Config& i_cfg) : cfg(i_cfg), iohandler(cfg),
                                            tsmanager(iohandler, cfg), gridsmanager(iohandler, cfg), interpolator(cfg, tsmanager, gridsmanager),
                                            downscaling(false), virtual_stations(false)
{
	initIOManager();
}

void IOManager::initIOManager()
{
	cfg.getValue("Virtual_stations", "Input", virtual_stations, IOUtils::nothrow);
	cfg.getValue("Downscaling", "Input", downscaling, IOUtils::nothrow);
}

void IOManager::setProcessingLevel(const unsigned int& i_level)
{
	tsmanager.setProcessingLevel(i_level);
	gridsmanager.setProcessingLevel(i_level);
}

void IOManager::setMinBufferRequirements(const double& buffer_size, const double& buff_before)
{
	tsmanager.setMinBufferRequirements(buffer_size, buff_before);
}

double IOManager::getAvgSamplingRate() const
{
	return tsmanager.getAvgSamplingRate();
}

const Config IOManager::getConfig() const
{
	return cfg;
}

void IOManager::push_meteo_data(const IOUtils::ProcessingLevel& level, const Date& date_start, const Date& date_end,
                                const std::vector< METEO_SET >& vecMeteo)
{
	tsmanager.push_meteo_data(level, date_start, date_end, vecMeteo);
}

size_t IOManager::getStationData(const Date& date, STATIONS_SET& vecStation)
{
	return tsmanager.getStationData(date, vecStation);
}

//for an interval of data: decide whether data should be filtered or raw
size_t IOManager::getMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< METEO_SET >& vecVecMeteo)
{
	return tsmanager.getMeteoData(dateStart, dateEnd, vecVecMeteo); //equivalent with the number of stations that have data
}

void IOManager::clear_cache()
{
	tsmanager.clear_cache();
	gridsmanager.clear_cache();
}

void IOManager::add_to_points_cache(const Date& i_date, const METEO_SET& vecMeteo)
{
	tsmanager.add_to_points_cache(i_date, vecMeteo);
}

//data can be raw or processed (filtered, resampled)
size_t IOManager::getMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();

	if(!virtual_stations && !downscaling) { //this is the usual case
		tsmanager.getMeteoData(i_date, vecMeteo);
	} else {
		if (virtual_stations)
			interpolator.getVirtualMeteoData(Meteo2DInterpolator::VSTATIONS, i_date, vecMeteo);
		if (downscaling)
			interpolator.getVirtualMeteoData(Meteo2DInterpolator::SMART_DOWNSCALING, i_date, vecMeteo);
	}

	return vecMeteo.size();
}

void IOManager::writeMeteoData(const std::vector< METEO_SET >& vecMeteo, const std::string& name)
{
	tsmanager.writeMeteoData(vecMeteo, name);
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result)
{
	string info_string;
	const bool status = getMeteoData(date, dem, meteoparam, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
	return status;
}

bool IOManager::getMeteoData(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                  Grid2DObject& result, std::string& info_string)
{
	interpolator.interpolate(date, dem, meteoparam, result, info_string);
	return (!result.empty());
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result)
{
	string info_string;
	interpolate(date, dem, meteoparam, in_coords, result, info_string);
	cerr << "[i] Interpolating " << MeteoData::getParameterName(meteoparam);
	cerr << " (" << info_string << ") " << endl;
}

void IOManager::interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string)
{
	interpolator.interpolate(date, dem, meteoparam, in_coords, result, info_string);
}

void IOManager::read2DGrid(Grid2DObject& grid2D, const std::string& filename)
{
	gridsmanager.read2DGrid(grid2D, filename);
}

void IOManager::read2DGrid(Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	gridsmanager.read2DGrid(grid2D, parameter, date);
}

void IOManager::readDEM(DEMObject& grid2D)
{
	gridsmanager.readDEM(grid2D);
}

void IOManager::readLanduse(Grid2DObject& grid2D)
{
	gridsmanager.readLanduse(grid2D);
}

void IOManager::readAssimilationData(const Date& date, Grid2DObject& grid2D)
{
	gridsmanager.readAssimilationData(date, grid2D);
}

void IOManager::readPOI(std::vector<Coords>& cpa)
{
	iohandler.readPOI(cpa);
}

void IOManager::write2DGrid(const Grid2DObject& grid2D, const std::string& name)
{
	gridsmanager.write2DGrid(grid2D, name);
}

void IOManager::write2DGrid(const Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	gridsmanager.write2DGrid(grid2D, parameter, date);
}

const std::string IOManager::toString() const {
	ostringstream os;
	os << "<IOManager>\n";
	os << "Config cfg = " << hex << &cfg << dec << "\n";
	os << iohandler.toString();
	os << tsmanager.toString();
	os << gridsmanager.toString();
	os << interpolator.toString();
	os << "Downscaling = " << downscaling << "\n";
	os << "Virtual stations = " << virtual_stations << "\n";
	os << "</IOManager>\n";
	return os.str();
}

} //namespace
