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
#include <meteoio/MeteoProcessor.h>

#include <algorithm>

using namespace std;

namespace mio {

MeteoProcessor::MeteoProcessor(const Config& cfg) : mi1d(cfg), processing_stack()
{
	//Parse [Filters] section, create processing stack for each configured parameter
	set<string> set_of_used_parameters;
	getParameters(cfg, set_of_used_parameters);

	for (set<string>::const_iterator it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it){
		ProcessingStack* tmp = new ProcessingStack(cfg, *it);
		processing_stack[*it] = tmp;
	}
}

MeteoProcessor::~MeteoProcessor()
{
	//clean up heap memory
	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); ++it)
		delete it->second;
}

void MeteoProcessor::getParameters(const Config& cfg, std::set<std::string>& set_parameters)
{
	std::vector<std::string> vec_keys;
	cfg.findKeys(vec_keys, std::string(), "Filters");

	for (size_t ii=0; ii<vec_keys.size(); ++ii){
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			const string tmp = vec_keys[ii].substr(0,found);
			set_parameters.insert(tmp);
		}
	}
}

void MeteoProcessor::getWindowSize(ProcessingProperties& o_properties) const
{
	ProcessingProperties tmp;

	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); ++it){
		(*(it->second)).getWindowSize(tmp);
		compareProperties(tmp, o_properties);
	}

	//Also take the Meteo1DInterpolator into account:
	mi1d.getWindowSize(tmp);
	compareProperties(tmp, o_properties);
}

void MeteoProcessor::compareProperties(const ProcessingProperties& newprop, ProcessingProperties& current)
{
	current.points_before = max(current.points_before, newprop.points_before);
	current.points_after = max(current.points_after, newprop.points_after);

	if (newprop.time_before > current.time_before)
		current.time_before = newprop.time_before;

	if (newprop.time_after > current.time_after)
		current.time_after = newprop.time_after;
}

void MeteoProcessor::process(const std::vector< std::vector<MeteoData> >& ivec,
                             std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	//call the different processing stacks
	std::vector< std::vector<MeteoData> > vec_tmp;

	for (map<string, ProcessingStack*>::const_iterator it=processing_stack.begin(); it != processing_stack.end(); ++it){
		if (it==processing_stack.begin()){
			(*(it->second)).process(ivec, ovec, second_pass);
		} else {
			(*(it->second)).process(vec_tmp, ovec, second_pass);
		}
		vec_tmp = ovec;
	}

	if (processing_stack.empty())
		ovec = ivec;
}

bool MeteoProcessor::resample(const Date& date, const std::vector<MeteoData>& ivec, MeteoData& md)
{
	return mi1d.resampleData(date, ivec, md);
}

const std::string MeteoProcessor::toString() const {
	std::ostringstream os;
	os << "<MeteoProcessor>\n";
	os << mi1d.toString();
	os << "Processing stacks:\n";
	map<string, ProcessingStack*>::const_iterator it;
	for (it=processing_stack.begin(); it != processing_stack.end(); ++it){
		//os << setw(10) << it->first.toString() << "::"; //the processing stack already contains it
		os << (*it->second).toString();
	}
	os << "</MeteoProcessor>\n";
	return os.str();
}

} //namespace
