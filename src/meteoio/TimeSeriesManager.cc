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

#include <meteoio/TimeSeriesManager.h>

#include <algorithm>

using namespace std;

namespace mio {

TimeSeriesManager::TimeSeriesManager(IOHandler& in_iohandler, const Config& in_cfg) : cfg(in_cfg), iohandler(in_iohandler),
                                            meteoprocessor(in_cfg), dataGenerator(in_cfg),
                                            proc_properties(), point_cache(), raw_buffer(), filtered_cache(),
                                            chunk_size(), buff_before(),
                                            processing_level(IOUtils::filtered | IOUtils::resampled | IOUtils::generated)
{
	meteoprocessor.getWindowSize(proc_properties);
	setDfltBufferProperties();
}

void TimeSeriesManager::setDfltBufferProperties()
{
	double chunk_size_days = 15.; //default chunk size value
	cfg.getValue("BUFF_CHUNK_SIZE", "General", chunk_size_days, IOUtils::nothrow); //in days
	chunk_size = Duration(chunk_size_days, 0);

	//get buffer centering options
	double buff_centering = -1.;
	double buff_start = -1.;
	cfg.getValue("BUFF_CENTERING", "General", buff_centering, IOUtils::nothrow);
	cfg.getValue("BUFF_BEFORE", "General", buff_start, IOUtils::nothrow);
	if ((buff_centering != -1.) && (buff_start != -1.))
		throw InvalidArgumentException("Please do NOT provide both BUFF_CENTERING and BUFF_BEFORE!!", AT);

	if (buff_start != -1.){
		buff_before = Duration(buff_start, 0);
	} else {
		if (buff_centering != -1.){
			if ((buff_centering < 0.) || (buff_centering > 1.))
				throw InvalidArgumentException("BUFF_CENTERING must be between 0 and 1", AT);

			buff_before = chunk_size * buff_centering;
		} else {
			buff_before = chunk_size * 0.1; //10% centering by default
		}
	}

	//if buff_before>chunk_size, we will have a problem (ie: we won't ever read the whole data we need)
	if(buff_before>chunk_size) chunk_size = buff_before;
	//BUG: if we do this, we still have the meteo1d window in the way
	//-> we end up not reading enough data and rebuffering...
}

void TimeSeriesManager::setMinBufferRequirements(const double& i_chunk_size, const double& i_buff_before)
{
	if(i_buff_before!=IOUtils::nodata) {
		const Duration app_buff_before(i_buff_before, 0);
		if(app_buff_before>buff_before) buff_before = app_buff_before;
	}
	if(i_chunk_size!=IOUtils::nodata) {
		const Duration app_chunk_size(i_chunk_size, 0);
		if(app_chunk_size>chunk_size) chunk_size = app_chunk_size;
	}

	//if buff_before>chunk_size, we will have a problem (ie: we won't ever read the whole data we need)
	if(buff_before>chunk_size) chunk_size = buff_before;
}

void TimeSeriesManager::setProcessingLevel(const unsigned int& i_level)
{
	if (i_level >= IOUtils::num_of_levels)
		throw InvalidArgumentException("The processing level is invalid", AT);

	if (((i_level & IOUtils::raw) == IOUtils::raw)
	    && ((i_level & IOUtils::filtered) == IOUtils::filtered))
		throw InvalidArgumentException("The processing level is invalid (raw and filtered at the same time)", AT);

	processing_level = i_level;
}

void TimeSeriesManager::push_meteo_data(const IOUtils::ProcessingLevel& level, const Date& date_start, const Date& date_end,
                                const std::vector< METEO_SET >& vecMeteo)
{
	//perform check on date_start and date_end
	if (date_end < date_start) {
		std::ostringstream ss;
		ss << "Trying to push data set from " << date_start.toString(Date::ISO) << " to " << date_end.toString(Date::ISO) << ". ";
		ss << " Obviously, date_start should be less than date_end!";
		throw InvalidArgumentException(ss.str(), AT);
	}

	if (level == IOUtils::filtered) {
		filtered_cache.push(date_start, date_end, vecMeteo);
	} else if (level == IOUtils::raw) {
		filtered_cache.clear();
		raw_buffer.push(date_start, date_end, vecMeteo);
	} else {
		throw InvalidArgumentException("The processing level is invalid (should be raw OR filtered)", AT);
	}

	point_cache.clear(); //clear point cache, so that we don't return resampled values of deprecated data
}

size_t TimeSeriesManager::getStationData(const Date& date, STATIONS_SET& vecStation)
{
	vecStation.clear();

	if (processing_level == IOUtils::raw){
		iohandler.readStationData(date, vecStation);
	} else {
		iohandler.readStationData(date, vecStation);
	}

	return vecStation.size();
}

//for an interval of data: decide whether data should be filtered or raw
size_t TimeSeriesManager::getMeteoData(const Date& dateStart, const Date& dateEnd, std::vector< METEO_SET >& vecVecMeteo)
{
	vecVecMeteo.clear();
	if (processing_level == IOUtils::raw){
		iohandler.readMeteoData(dateStart, dateEnd, vecVecMeteo);
	} else {
		const bool success = filtered_cache.get(dateStart, dateEnd, vecVecMeteo);

		if (!success){
			vector< vector<MeteoData> > tmp_meteo;
			fillRawBuffer(dateStart, dateEnd);
			raw_buffer.get(dateStart, dateEnd, tmp_meteo);

			//now it needs to be secured that the data is actually filtered, if configured
			if ((IOUtils::filtered & processing_level) == IOUtils::filtered){
				fill_filtered_cache();
				filtered_cache.get(dateStart, dateEnd, vecVecMeteo);
			} else {
				vecVecMeteo = tmp_meteo;
			}
		}

		if ((IOUtils::generated & processing_level) == IOUtils::generated){
			dataGenerator.createParameters(vecVecMeteo);
			dataGenerator.fillMissing(vecVecMeteo);
		}
	}

	return vecVecMeteo.size(); //equivalent with the number of stations that have data
}

size_t TimeSeriesManager::getMeteoData(const Date& i_date, METEO_SET& vecMeteo)
{
	vecMeteo.clear();

	//1. Check whether user wants raw data or processed data
	//The first case: we are looking at raw data directly, only unresampled values are considered, exact date match
	if (processing_level == IOUtils::raw) {
		vector< vector<MeteoData> > vec_cache;
		const Duration eps(1./(24.*3600.), 0.);
		iohandler.readMeteoData(i_date-eps, i_date+eps, vec_cache);
		for (size_t ii=0; ii<vec_cache.size(); ii++){ //for every station
			const size_t index = IOUtils::seek(i_date, vec_cache[ii], true);
			if (index != IOUtils::npos)
				vecMeteo.push_back(vec_cache[ii][index]); //Insert station into vecMeteo
		}
		return vecMeteo.size();
	}

	//2.  Check which data point is available, buffered locally
	const map<Date, vector<MeteoData> >::const_iterator it = point_cache.find(i_date);
	if (it != point_cache.end()){
		vecMeteo = it->second;
		return vecMeteo.size();
	}

	//Let's make sure we have the data we need, in the filtered_cache or in vec_cache
	const Date buffer_start( i_date-proc_properties.time_before ), buffer_end( i_date+proc_properties.time_after );
	vector< vector<MeteoData> >* data = NULL; //reference to either filtered_cache or raw_buffer
	if ((IOUtils::filtered & processing_level) == IOUtils::filtered){
		const bool cached = (!filtered_cache.empty()) && (filtered_cache.getBufferStart() <= buffer_start) && (filtered_cache.getBufferEnd() >= buffer_end);
		if (!cached) {
			//explicit caching, rebuffer if necessary
			fillRawBuffer(buffer_start, buffer_end);
			fill_filtered_cache();
		}
		data = &filtered_cache.getBuffer();
	} else { //data to be resampled should be IOUtils::raw
		fillRawBuffer(buffer_start, buffer_end);
		data = &raw_buffer.getBuffer();
	}

	if ((IOUtils::resampled & processing_level) != IOUtils::resampled) { //no resampling required
		for (size_t ii=0; ii<(*data).size(); ii++) { //for every station
			const size_t index = IOUtils::seek(i_date, (*data)[ii], true); //needs to be an exact match
			if (index != IOUtils::npos)
				vecMeteo.push_back((*data)[ii][index]); //Insert station into vecMeteo
		}
	} else { //resampling required
		MeteoData md;
		for (size_t ii=0; ii<(*data).size(); ii++) { //for every station
			const bool success = meteoprocessor.resample(i_date, (*data)[ii], md);
			if (success) vecMeteo.push_back(md);
		}
	}

	if ((IOUtils::generated & processing_level) == IOUtils::generated) {
		dataGenerator.createParameters(vecMeteo);
		dataGenerator.fillMissing(vecMeteo);
	}

	add_to_points_cache(i_date, vecMeteo); //Store result in the local cache

	return vecMeteo.size();
}

void TimeSeriesManager::writeMeteoData(const std::vector< METEO_SET >& vecMeteo, const std::string& name)
{
	if (processing_level == IOUtils::raw){
		iohandler.writeMeteoData(vecMeteo, name);
	} else {
		iohandler.writeMeteoData(vecMeteo, name);
	}
}

double TimeSeriesManager::getAvgSamplingRate() const
{
	return raw_buffer.getAvgSamplingRate();
}

/**
 * @brief Filter the whole raw meteo data buffer
 */
void TimeSeriesManager::fill_filtered_cache()
{
	if ((IOUtils::filtered & processing_level) == IOUtils::filtered){
		filtered_cache.clear(); //HACK until we get true ringbuffers, to prevent eating up all memory
		meteoprocessor.process(raw_buffer.getBuffer(), filtered_cache.getBuffer());
		filtered_cache.setBufferStart(raw_buffer.getBufferStart()); //HACK
		filtered_cache.setBufferEnd(raw_buffer.getBufferEnd()); //HACK
	}
}

void TimeSeriesManager::add_to_points_cache(const Date& i_date, const METEO_SET& vecMeteo)
{
	//Check cache size, delete oldest elements if necessary
	if (point_cache.size() > 2000) {
		point_cache.clear(); //HACK: implement a true ring buffer!
	}

	point_cache[i_date] = vecMeteo;
}

void TimeSeriesManager::clear_cache()
{
	raw_buffer.clear();
	filtered_cache.clear();
	point_cache.clear();
}

void TimeSeriesManager::fillRawBuffer(const Date& date_start, const Date& date_end)
{
	const Date new_start( date_start-buff_before ); //taking centering into account
	const Date new_end( max(new_start + chunk_size, date_end) );

	raw_buffer.clear(); //HACK until we have a proper ring buffer to avoid eating up all memory...

	if (raw_buffer.empty()) {
		vector< METEO_SET > vecMeteo;
		iohandler.readMeteoData(new_start, new_end, vecMeteo);
		raw_buffer.push(new_start, new_end, vecMeteo);
		return;
	}

	const Date buffer_start( raw_buffer.getBufferStart() );
	const Date buffer_end( raw_buffer.getBufferEnd() );
	if (new_start>buffer_end || new_end<buffer_start) { //easy: full rebuffer
		vector< METEO_SET > vecMeteo;
		iohandler.readMeteoData(new_start, new_end, vecMeteo);
		raw_buffer.push(new_start, new_end, vecMeteo);
		return;
	}

	if (new_start<buffer_start) { //some data must be inserted before
		vector< METEO_SET > vecMeteo;
		iohandler.readMeteoData(new_start, buffer_start, vecMeteo);
		raw_buffer.push(new_start, buffer_start, vecMeteo);
	}

	if (new_end>buffer_end) { //some data must be inserted after. Keep in mind both before and after could happen simultaneously!
		vector< METEO_SET > vecMeteo;
		iohandler.readMeteoData(buffer_end, new_end, vecMeteo);
		raw_buffer.push(buffer_end, new_end, vecMeteo);
	}

}

const std::string TimeSeriesManager::toString() const {
	ostringstream os;
	os << "<TimeSeriesManager>\n";
	os << "Config& cfg = " << hex << &cfg << dec << "\n";
	os << "IOHandler& iohandler = " << hex << &iohandler << dec << "\n";
	os << meteoprocessor.toString();
	os << "Processing level = " << processing_level << "\n";
	os << dataGenerator.toString();

	os << "RawBuffer:\n" << raw_buffer.toString();
	os << "Filteredcache:\n" << raw_buffer.toString();

	//display point_cache
	size_t count=0;
	size_t min_stations=std::numeric_limits<size_t>::max();
	size_t max_stations=0;
	std::map<Date, std::vector<MeteoData> >::const_iterator iter = point_cache.begin();
	for (; iter != point_cache.end(); ++iter) {
		const size_t nb_stations = iter->second.size();
		if(nb_stations>max_stations) max_stations=nb_stations;
		if(nb_stations<min_stations) min_stations=nb_stations;
		count++;
	}

	if(count==0) {
		os << "Resampled cache is empty\n";
	}
	if(count==1) {
		os << "Resampled cache content (";
		if(max_stations==min_stations)
			os << min_stations;
		else
			os << min_stations << " to " << max_stations;
		os << " station(s))\n";
		os << std::setw(22) << point_cache.begin()->first.toString(Date::ISO) << " - 1 timestep\n";
	}
	if(count>1) {
		const double avg_sampling = ( (point_cache.rbegin()->first.getJulian()) - (point_cache.begin()->first.getJulian()) ) / (double)(count-1);

		os << "Resampled cache content (";
		if(max_stations==min_stations)
			os << min_stations;
		else
			os << min_stations << " to " << max_stations;
		os << " station(s))\n";
		os << std::setw(22) << point_cache.begin()->first.toString(Date::ISO);
		os << " - " << point_cache.rbegin()->first.toString(Date::ISO);
		os << " - " << count << " timesteps (" << setprecision(3) << fixed << avg_sampling*24.*3600. << " s sampling rate)";
	}

	os << "</TimeSeriesManager>\n";
	return os.str();
}

} //namespace
