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

#include <meteoio/dataClasses/Buffer.h>

#include <algorithm>

using namespace std;

namespace mio {

bool MeteoBuffer::get(const Date& date, METEO_SET &vecMeteo) const
{
	vecMeteo.clear();

	if (empty() || (date < ts_start) || (date > ts_end)) //data is NOT fully in cache
		return false;

	for (size_t ii=0; ii<ts_buffer.size(); ii++) { //loop over stations
		if (ts_buffer[ii].empty()) continue; //no data in buffer for this station

		const size_t pos = IOUtils::seek(date, ts_buffer[ii], true);
		if (pos!=IOUtils::npos)
			vecMeteo.push_back( ts_buffer[ii][pos] );
	}

	return true;
}

bool MeteoBuffer::get(const Date& date_start, const Date& date_end, std::vector< METEO_SET > &vecMeteo) const
{
	vecMeteo.clear();

	if (empty() || (date_start < ts_start) || (date_end > ts_end)) //data is NOT fully in cache
		return false;

	for (size_t ii=0; ii<ts_buffer.size(); ii++) { //loop over stations
		vecMeteo.push_back( vector<MeteoData>() );  //insert one empty vector of MeteoData

		if (ts_buffer[ii].empty()) continue; //no data in buffer for this station
		if (ts_buffer[ii].front().date>date_end || ts_buffer[ii].back().date<date_start) continue; //no data in buffer for this station

		size_t pos_start = IOUtils::seek(date_start, ts_buffer[ii], false);
		if (pos_start==IOUtils::npos) pos_start = 0;
		size_t pos_end = IOUtils::seek(date_end, ts_buffer[ii], false);
		if (pos_end==IOUtils::npos) pos_end = ts_buffer[ii].size() - 1;

		vecMeteo[ii].reserve(pos_end-pos_start+1); //weird that the "insert" does not handle it internally...
		vecMeteo[ii].insert(vecMeteo[ii].begin(), ts_buffer[ii].begin()+pos_start, ts_buffer[ii].begin()+pos_end);
	}

	return true;
}

bool MeteoBuffer::empty() const
{ //the ts_buffer could be empty if there was no data between the provided dates
  //so the empty criteria is if the ts_start and ts_end are valid or undef
	return (ts_start.isUndef() || ts_end.isUndef());
}

void MeteoBuffer::clear()
{
	ts_buffer.clear();
	ts_start.setUndef(true);
	ts_end.setUndef(true);
}

void MeteoBuffer::push(const Date& date_start, const Date& date_end, const std::vector< METEO_SET >& vecMeteo)
{
	if (empty()) {
		ts_start = date_start;
		ts_end = date_end;
		ts_buffer = vecMeteo;
		return;
	}

	const size_t nrStationsBuffer = ts_buffer.size();
	const size_t nrStationsPush = vecMeteo.size();

	//check that we are dealing with the same stations
	if (nrStationsBuffer!=nrStationsPush) {
		ostringstream ss;
		ss << "The number of stations changed over time from " << nrStationsBuffer << " to " << nrStationsPush << ", ";
		ss << "this is not handled yet!";
		throw IOException(ss.str(), AT);
	}

	for (size_t ii=0; ii<nrStationsBuffer; ii++) { //for all stations
		if (ts_buffer[ii].empty() || vecMeteo[ii].empty())
			continue;
		if (ts_buffer[ii].front().meta.getHash()!=vecMeteo[ii].front().meta.getHash()) {
			ostringstream ss;
			ss << "The stations changed over time from " << ts_buffer[ii].front().meta.getHash() << " to " << vecMeteo[ii].front().meta.getHash() << ", ";
			ss << "this is not handled yet!";
			throw IOException(ss.str(), AT);
		}
	}

	//now, do the append/merge
	for (size_t ii=0; ii<nrStationsBuffer; ii++) { //for all stations
		if (vecMeteo[ii].empty()) continue;

		if (ts_buffer[ii].empty()) {
			ts_buffer[ii] = vecMeteo[ii];
			continue;
		}

		const Date buffer_start = ts_buffer[ii].front().date;
		const Date buffer_end = ts_buffer[ii].back().date;
		const Date data_start = vecMeteo[ii].front().date;
		const Date data_end = vecMeteo[ii].back().date;

		if (data_start>buffer_end) { //the data simply fits to the end
			ts_buffer[ii].insert(ts_buffer[ii].end(), vecMeteo[ii].begin(), vecMeteo[ii].begin()+vecMeteo[ii].size());
			continue;
		}
		if (data_end<buffer_start) { //the data simply fits at the start
			ts_buffer[ii].insert(ts_buffer[ii].begin(), vecMeteo[ii].begin(), vecMeteo[ii].begin()+vecMeteo[ii].size());
			continue;
		}

		//there is some overlap, only copy data that does NOT overlap
		if (data_start<buffer_start) {
			const size_t pos = IOUtils::seek(buffer_start, vecMeteo[ii], false); //returns the first date >=
			ts_buffer[ii].insert(ts_buffer[ii].begin(), vecMeteo[ii].begin(), vecMeteo[ii].begin()+pos-1);
		}
		if (data_end>buffer_end) {
			size_t pos = IOUtils::seek(buffer_end, vecMeteo[ii], false); //returns the first date >=
			if (vecMeteo[ii][pos].date == buffer_end) pos++; //to make sure we have an element that is not already in buffer
			ts_buffer[ii].insert(ts_buffer[ii].end(), vecMeteo[ii].begin()+pos, vecMeteo[ii].begin()+vecMeteo[ii].size());
		}
	}

	ts_start = min(ts_start, date_start);
	ts_end = max(ts_end, date_end);
}

/*void MeteoBuffer::push(const Date& date, const METEO_SET& vecMeteo)
{

}*/

double MeteoBuffer::getAvgSamplingRate() const
{
	if (ts_buffer.empty())
		return IOUtils::nodata;

	const size_t nr_stations = ts_buffer.size();
	double sum = 0;
	for (size_t ii=0; ii<nr_stations; ii++){ //loop over all stations
		if(!ts_buffer[ii].empty()) {
			const std::vector<MeteoData>& curr_station = ts_buffer[ii];
			const double days = curr_station.back().date.getJulian() - curr_station.front().date.getJulian();

			//add the average sampling rate for this station
			const size_t nr_data_pts = ts_buffer[ii].size();
			if(days>0.) sum += (double)(nr_data_pts-1) / days; //the interval story: 2 points define 1 interval!
		}
	}
	if (sum > 0.){
		return ((double)sum / (double)(nr_stations*24*3600)); //in points per seconds, ie Hz
	}

	return IOUtils::nodata;
}

Date MeteoBuffer::getBufferStart() const
{
	return ts_start;
}

Date MeteoBuffer::getBufferEnd() const
{
	return ts_end;
}

std::vector< METEO_SET >& MeteoBuffer::getBuffer()
{
	return ts_buffer;
}

void MeteoBuffer::setBufferStart(const Date& date) {
	ts_start = date;
}

void MeteoBuffer::setBufferEnd(const Date& date) {
	ts_end = date;
}

const std::string MeteoBuffer::toString() const
{
	ostringstream os;
	os << "<MeteoBuffer>\n";

	os << "Buffer content (" << ts_buffer.size() << " stations)\n";
	for(size_t ii=0; ii<ts_buffer.size(); ii++) {
		if (!ts_buffer[ii].empty()){
			os << std::setw(10) << ts_buffer[ii].front().meta.stationID << " = "
			   << ts_buffer[ii].front().date.toString(Date::ISO) << " - "
			   << ts_buffer[ii].back().date.toString(Date::ISO) << ", "
			   << ts_buffer[ii].size() << " timesteps" << endl;
		}
	}

	os << "</MeteoBuffer>\n";
	return os.str();
}

/********************************************************************************************/
/****************************** GridBuffer class ***********************************************/
/********************************************************************************************/

GridBuffer::GridBuffer(const size_t& in_max_grids) 
           : mapBufferedGrids(), mapBufferedDEMs(), mapBufferedInfos(), IndexBufferedGrids(), IndexBufferedDEMs(), max_grids(in_max_grids) 
{}


bool GridBuffer::get(Grid2DObject& grid, const std::string& grid_hash) const
{
	if (IndexBufferedGrids.empty())
		return false;

	const std::map<std::string, Grid2DObject>::const_iterator it = mapBufferedGrids.find( grid_hash );
	if (it != mapBufferedGrids.end()) { //already in map
		grid = (*it).second;
		return true;
	}

	return false;
}

bool GridBuffer::get(Grid2DObject& grid, const std::string& grid_hash, std::string& grid_info) const
{
	grid_info.clear();
	if (max_grids==0) return false;
	
	if (get(grid, grid_hash)==true) {
		const std::map<std::string, std::string>::const_iterator it = mapBufferedInfos.find( grid_hash );
		if (it != mapBufferedInfos.end()) { //found in map
			grid_info = (*it).second;
		}
		return true;
	}
	
	return false;
}

bool GridBuffer::get(Grid2DObject& grid, const MeteoGrids::Parameters& parameter, const Date& date) const
{
	const string grid_hash = date.toString(Date::ISO)+"::"+MeteoGrids::getParameterName(parameter);
	return get(grid, grid_hash);
}

bool GridBuffer::get(DEMObject& grid, const std::string& grid_hash) const
{
	if (IndexBufferedDEMs.empty())
		return false;

	const std::map<std::string, DEMObject>::const_iterator it = mapBufferedDEMs.find( grid_hash );
	if (it != mapBufferedDEMs.end()) { //already in map
		//properties of the passed dem
		const DEMObject::update_type in_ppt = (DEMObject::update_type)grid.getUpdatePpt();
		const DEMObject::slope_type in_slope_alg = (DEMObject::slope_type)grid.getDefaultAlgorithm();

		//get cached dem and its properties
		grid = (*it).second;
		const DEMObject::update_type buff_ppt = (DEMObject::update_type)grid.getUpdatePpt();
		const DEMObject::slope_type buff_slope_alg = (DEMObject::slope_type)grid.getDefaultAlgorithm();

		if (in_ppt!=buff_ppt || in_slope_alg!=buff_slope_alg) {
			//if passed properties != cached properties -> set the returned one to the passed properties
			grid.setDefaultAlgorithm(in_slope_alg);
			grid.setUpdatePpt(in_ppt);
			grid.update();
		}

		return true;
	}

	return false;
}

void GridBuffer::push(const DEMObject& grid, const std::string& grid_hash)
{
	if (max_grids==0) return;

	if(IndexBufferedDEMs.size() >= max_grids) { //we need to remove the oldest grid
		const string tmp_hash = IndexBufferedDEMs.front();
		mapBufferedDEMs.erase( mapBufferedDEMs.find( tmp_hash ) );
		//swap followed by pop_back is faster than erase()
		swap( IndexBufferedDEMs.front(), IndexBufferedDEMs.back() );
		IndexBufferedDEMs.pop_back();
	}
	mapBufferedDEMs[ grid_hash ] = grid;
	IndexBufferedDEMs.push_back( grid_hash );
}

void GridBuffer::push(const Grid2DObject& grid, const std::string& grid_hash, const std::string& grid_info)
{
	if (max_grids==0) return;

	if(IndexBufferedGrids.size() >= max_grids) { //we need to remove the oldest grid
		const string tmp_hash = IndexBufferedGrids.front();
		mapBufferedGrids.erase( mapBufferedGrids.find( tmp_hash ) );
		mapBufferedInfos.erase( mapBufferedInfos.find( tmp_hash ) );
		//swap followed by pop_back is faster than erase()
		swap( IndexBufferedGrids.front(), IndexBufferedGrids.back() );
		IndexBufferedGrids.pop_back();
	}
	mapBufferedGrids[ grid_hash ] = grid;
	mapBufferedInfos[ grid_hash ] = grid_info;
	IndexBufferedGrids.push_back( grid_hash );
}

void GridBuffer::push(const Grid2DObject& grid, const std::string& grid_hash)
{
	if (max_grids==0) return;
	
	push(grid, grid_hash, "");
}

void GridBuffer::push(const Grid2DObject& grid, const MeteoGrids::Parameters& parameter, const Date& date)
{
	if (max_grids==0) return;
	
	const string grid_hash = date.toString(Date::ISO)+"::"+MeteoGrids::getParameterName(parameter);
	push(grid, grid_hash, "");
}

const std::string GridBuffer::toString() const
{
	ostringstream os;
	os << "<GridBuffer>\n";
	os << "Max buffered grids = " << max_grids << "\n";

	//cache content
	os << "Cached grids: " << mapBufferedGrids.size() << "\n";
	std::map<std::string, Grid2DObject>::const_iterator it_grid;
	for (it_grid=mapBufferedGrids.begin(); it_grid != mapBufferedGrids.end(); ++it_grid){
		os << setw(10) << "Grid " << it_grid->first << "\n";
	}
	
	//dem buffer
	os << "Cached dems: " << mapBufferedDEMs.size() << "\n";
	std::map<std::string, DEMObject>::const_iterator it_dem;
	for (it_dem=mapBufferedDEMs.begin(); it_dem != mapBufferedDEMs.end(); ++it_dem){
		os << setw(10) << "Dem " << it_dem->first << "\n";
	}

	os << "</GridBuffer>\n";
	return os.str();
}


} //end namespace
