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
#include <meteoio/Meteo1DInterpolator.h>

using namespace std;

namespace mio {

Meteo1DInterpolator::Meteo1DInterpolator(const Config& in_cfg)
                     : cfg(in_cfg), window_size(86400.),
                       mapAlgorithms()
{
	//default window_size is 2 julian days
	cfg.getValue("WINDOW_SIZE", "Interpolations1D", window_size, IOUtils::nothrow);
	if (window_size <= 1.)
		throw IOException("WINDOW_SIZE not valid, it should be a duration in seconds at least greater than 1", AT);
	window_size /= 86400.; //user uses seconds, internally julian day is used

	//read the Config object to create the resampling algorithms for each
	//MeteoData::Parameters parameter (i.e. each member variable like ta, p, psum, ...)
	for (size_t ii=0; ii<MeteoData::nrOfParameters; ii++){ //loop over all MeteoData member variables
		const std::string parname = MeteoData::getParameterName(ii); //Current parameter name
		vector<string> vecArgs;
		const string algo_name = getInterpolationForParameter(parname, vecArgs);
		mapAlgorithms[parname] = ResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname, window_size, vecArgs);
	}
}

Meteo1DInterpolator::~Meteo1DInterpolator()
{
	map< string, ResamplingAlgorithms* >::iterator it;
	for (it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it)
		delete it->second;
}

void Meteo1DInterpolator::getWindowSize(ProcessingProperties& o_properties) const
{
	o_properties.points_before = 1;
	o_properties.points_after  = 1;
	o_properties.time_before   = Duration(window_size, 0.); //we will need to cut a window 2x larger so we can interpolate each point in the window
	o_properties.time_after    = Duration(window_size, 0.);
}

bool Meteo1DInterpolator::resampleData(const Date& date, const std::vector<MeteoData>& vecM, MeteoData& md)
{
	if (vecM.empty()) //Deal with case of the empty vector
		return false; //nothing left to do

	md = vecM.front(); //create a clone of one of the elements
	md.reset();   //set all values to IOUtils::nodata
	md.setDate(date);

	//Find element in the vector or the next index
	size_t index = IOUtils::seek(date, vecM, false);

	//Three cases
	ResamplingAlgorithms::ResamplingPosition elementpos = ResamplingAlgorithms::exact_match;
	if (index == IOUtils::npos) { //nothing found append new element at the left or right
		if (vecM.front().date > date) {
			elementpos = ResamplingAlgorithms::begin;
			index = 0;
		} else if (vecM.back().date < date) {
			elementpos = ResamplingAlgorithms::end;
			index = vecM.size() - 1;
		}
		md.setResampled(true);
	} else if ((index != IOUtils::npos) && (vecM[index].date != date)) {
		elementpos = ResamplingAlgorithms::before;
		md.setResampled(true);
	} else {
		md.setResampled(false);
	}

	//now, perform the resampling
	for (size_t ii=0; ii<md.getNrOfParameters(); ii++) {
		const std::string parname = md.getNameForParameter(ii); //Current parameter name
		const map< string, ResamplingAlgorithms* >::const_iterator it = mapAlgorithms.find(parname);
		if (it!=mapAlgorithms.end()) {
			it->second->resample(index, elementpos, ii, vecM, md);
		} else { //we are dealing with an extra parameter, we need to add it to the map first, so it will exist next time...
			vector<string> vecArgs;
			const string algo_name = getInterpolationForParameter(parname, vecArgs);
			mapAlgorithms[parname] = ResamplingAlgorithmsFactory::getAlgorithm(algo_name, parname, window_size, vecArgs);;
			mapAlgorithms[parname]->resample(index, elementpos, ii, vecM, md);
		}

		#ifdef DATA_QA
		const map< string, ResamplingAlgorithms* >::const_iterator it2 = mapAlgorithms.find(parname); //we have to re-find it in order to handle extra parameters
		if ((index != IOUtils::npos) && vecM[index](ii)!=md(ii)) {
			const string statName = md.meta.getStationName();
			const string statID = md.meta.getStationID();
			const string stat = (!statID.empty())? statID : statName;
			const string algo_name = it2->second->getAlgo();
			cout << "[DATA_QA] Resampling " << stat << "::" << parname << "::" << algo_name << " " << md.date.toString(Date::ISO_TZ) << " [" << md.date.toString(Date::ISO_WEEK) << "]\n";
		}
		#endif
	}

	return true; //successfull resampling
}

/**
 * @brief retrieve the resampling algorithm to be used for the 1D interpolation of meteo parameters.
 * The potential arguments are also extracted.
 * @param parname meteo parameter to deal with
 * @param vecArguments vector of arguments
 * @return algorithm name
 */
string Meteo1DInterpolator::getInterpolationForParameter(const std::string& parname, std::vector<std::string>& vecArguments) const
{
	string algo_name = "linear"; //the default resampling is linear
	cfg.getValue(parname+"::resample", "Interpolations1D", algo_name, IOUtils::nothrow);

	cfg.getValue(parname+"::"+algo_name, "Interpolations1D", vecArguments, IOUtils::nothrow);

	if (cfg.keyExists(parname+"::"+"args", "Interpolations1D")) //HACK: temporary until we consider everybody has migrated
		throw InvalidArgumentException("The syntax for Interpolations1D arguments has been changed. Please check the documentation and update your configuration file \""+cfg.getSourceName()+"\"", AT);

	return algo_name;
}

Meteo1DInterpolator& Meteo1DInterpolator::operator=(const Meteo1DInterpolator& source) {
	if (this != &source) {
		window_size = source.window_size;
		mapAlgorithms= source.mapAlgorithms;
	}
	return *this;
}

const std::string Meteo1DInterpolator::toString() const
{
	ostringstream os;
	os << "<Meteo1DInterpolator>\n";
	os << "Config& cfg = " << hex << &cfg << dec <<"\n";
	os << "Resampling algorithms:\n";
	map< string, ResamplingAlgorithms* >::const_iterator it;
	for (it=mapAlgorithms.begin(); it!=mapAlgorithms.end(); ++it) {
		//os << setw(10) << it->first << "::" << it->second->getAlgo() << "\n";
		os << it->second->toString() << "\n";
	}
	os << "</Meteo1DInterpolator>\n";

	return os.str();
}

} //namespace
