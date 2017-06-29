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

#include <meteoio/GridsManager.h>
#include <meteoio/dataClasses/Coords.h>

using namespace std;

namespace mio {

GridsManager::GridsManager(IOHandler& in_iohandler, const Config& in_cfg)
             : iohandler(in_iohandler), cfg(in_cfg), buffer(0),
               processing_level(IOUtils::filtered | IOUtils::resampled | IOUtils::generated)
{
	size_t max_grids = 10;
	cfg.getValue("BUFF_GRIDS", "General", max_grids, IOUtils::nothrow);
	buffer.setMaxGrids(max_grids);
}

/**
* @brief Set the desired ProcessingLevel
*        The processing level affects the way meteo data is read and processed
*        Three values are possible:
*        - IOUtils::raw data shall be read directly from the buffer
*        - IOUtils::filtered data shall be filtered before returned to the user
*        - IOUtils::resampled data shall be resampled before returned to the user
*          this only affects the function getMeteoData(const Date&, METEO_DATASET&);
*
*        The three values can be combined: e.g. IOUtils::filtered | IOUtils:resampled
* @param i_level The ProcessingLevel values that shall be used to process data
*/
void GridsManager::setProcessingLevel(const unsigned int& i_level)
{
	if (i_level >= IOUtils::num_of_levels)
		throw InvalidArgumentException("The processing level is invalid", AT);

	if (((i_level & IOUtils::raw) == IOUtils::raw)
	    && ((i_level & IOUtils::filtered) == IOUtils::filtered))
		throw InvalidArgumentException("The processing level is invalid (raw and filtered at the same time)", AT);

	processing_level = i_level;
}

void GridsManager::clear_cache()
{
	buffer.clear();
}

void GridsManager::read2DGrid(Grid2DObject& grid2D, const std::string& filename)
{
	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, filename);
	} else {
		if (buffer.get(grid2D, filename))
			return;

		iohandler.read2DGrid(grid2D, filename);
		buffer.push(grid2D, filename);
	}
}

void GridsManager::read2DGrid(Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (processing_level == IOUtils::raw){
		iohandler.read2DGrid(grid2D, parameter, date);
	} else {
		if (buffer.get(grid2D, parameter, date))
			return;

		iohandler.read2DGrid(grid2D, parameter, date);
		buffer.push(grid2D, parameter, date);
	}
}

//HACK buffer 3D grids!
void GridsManager::read3DGrid(Grid3DObject& grid_out, const std::string& i_filename)
{
	iohandler.read3DGrid(grid_out, i_filename);
}

void GridsManager::read3DGrid(Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.read3DGrid(grid_out, parameter, date);
}

void GridsManager::readDEM(DEMObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readDEM(grid2D);
	} else {
		if (buffer.get(grid2D, "/:DEM"))
			return;

		iohandler.readDEM(grid2D);
		buffer.push(grid2D, "/:DEM");
	}
}

void GridsManager::readLanduse(Grid2DObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readLanduse(grid2D);
	} else {
		if (buffer.get(grid2D, "/:LANDUSE"))
			return;

		iohandler.readLanduse(grid2D);
		buffer.push(grid2D, "/:LANDUSE");
	}
}

void GridsManager::readAssimilationData(const Date& date, Grid2DObject& grid2D)
{
	if (processing_level == IOUtils::raw){
		iohandler.readAssimilationData(date, grid2D);
	} else {
		const string grid_hash = "/:ASSIMILATIONDATA"+date.toString(Date::ISO);
		if (buffer.get(grid2D, grid_hash))
			return;

		iohandler.readAssimilationData(date, grid2D);
		buffer.push(grid2D, grid_hash);
	}
}

void GridsManager::write2DGrid(const Grid2DObject& grid2D, const std::string& name)
{
	iohandler.write2DGrid(grid2D, name);
}

void GridsManager::write2DGrid(const Grid2DObject& grid2D, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.write2DGrid(grid2D, parameter, date);
}

void GridsManager::write3DGrid(const Grid3DObject& grid_out, const std::string& options)
{
	iohandler.write3DGrid(grid_out, options);
}

void GridsManager::write3DGrid(const Grid3DObject& grid_out, const MeteoGrids::Parameters& parameter, const Date& date)
{
	iohandler.write3DGrid(grid_out, parameter, date);
}

const std::string GridsManager::toString() const {
	ostringstream os;
	os << "<GridsManager>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "IOHandler& iohandler = " << hex << &iohandler << dec << "\n";
	os << "Processing level = " << processing_level << "\n";
	os << buffer.toString();
	os << "</GridsManager>\n";
	return os.str();
}

} //namespace
