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
#include <meteoio/meteoFilters/ProcessingStack.h>

#include <algorithm>

using namespace std;

namespace mio {

ProcessingStack::ProcessingStack(const Config& cfg, const std::string& parname) : filter_stack(), param_name(parname)
{
	 //this is required by filters that need to read some parameters in extra files
	const string root_path = cfg.getConfigRootDir();
	vector<string> vecFilters;
	cfg.getValues(parname+"::filter", "Filters", vecFilters);

	const size_t nr_of_filters = vecFilters.size();
	for (size_t ii=0; ii<nr_of_filters; ii++){
		//create a processing block for each filter
		const string block_name = IOUtils::strToUpper( vecFilters[ii] );
		std::vector<std::string> vec_args;
		std::ostringstream tmp;
		tmp << param_name << "::arg" << (ii+1);

		//Read arguments
		cfg.getValue(tmp.str(), "Filters", vec_args, IOUtils::nothrow);
		filter_stack.push_back( BlockFactory::getBlock(block_name, vec_args, root_path) );
	}
}

ProcessingStack::~ProcessingStack()
{
	for (size_t ii=0; ii<filter_stack.size(); ii++)
		delete filter_stack[ii];
}

void ProcessingStack::getWindowSize(ProcessingProperties& o_properties)
{
	o_properties.points_before = 0;
	o_properties.points_after = 0;
	o_properties.time_after = Duration(0.0, 0.);
	o_properties.time_before = Duration(0.0, 0.);

	for (size_t jj=0; jj<filter_stack.size(); jj++){
		const ProcessingProperties& properties = (*filter_stack[jj]).getProperties();

		o_properties.points_before = max(o_properties.points_before, properties.points_before);
		o_properties.points_after = max(o_properties.points_after, properties.points_after);

		if (properties.time_before > o_properties.time_before)
			o_properties.time_before = properties.time_before;

		if (properties.time_after > o_properties.time_after)
			o_properties.time_after = properties.time_after;
	}
}

//this method applies the whole processing stack for all the stations, all the data points for one meteo param
//(as defined in the constructor)
void ProcessingStack::process(const std::vector< std::vector<MeteoData> >& ivec,
                              std::vector< std::vector<MeteoData> >& ovec, const bool& second_pass)
{
	const size_t nr_of_filters = filter_stack.size();
	const size_t nr_stations = ivec.size();
	ovec.resize( nr_stations );

	for (size_t ii=0; ii<nr_stations; ii++){ //for every station
		if ( ivec[ii].empty() ) continue; //no data, nothing to do!

		//pick one element and check whether the param_name parameter exists
		const size_t param = ivec[ii].front().getParameterIndex(param_name);
		if (param != IOUtils::npos){
			//since all filters start with ovec=ivec, maybe we could just swap pointers instead for copying
			std::vector<MeteoData> tmp( ivec[ii] );

			//Now call the filters one after another for the current station and parameter
			bool appliedFilter = false;
			for (size_t jj=0; jj<nr_of_filters; jj++){
				const ProcessingProperties::proc_stage& filter_stage = filter_stack[jj]->getProperties().stage;
				if ( second_pass && ((filter_stage==ProcessingProperties::first) || (filter_stage==ProcessingProperties::none)) )
					continue;

				if ( !second_pass && ((filter_stage==ProcessingProperties::second) || (filter_stage==ProcessingProperties::none)) )
					continue;

				appliedFilter = true;
				(*filter_stack[jj]).process(static_cast<unsigned int>(param), tmp, ovec[ii]);

				if (tmp.size() == ovec[ii].size()){
					#ifdef DATA_QA
					for (size_t kk=0; kk<ovec[ii].size(); kk++) {
						const double orig = tmp[kk](param);
						const double filtered = ovec[ii][kk](param);
						if (orig!=filtered) {
							const string statName = ovec[ii][kk].meta.getStationName();
							const string statID = ovec[ii][kk].meta.getStationID();
							const string stat = (!statID.empty())? statID : statName;
							const string filtername = (*filter_stack[jj]).getName();
							cout << "[DATA_QA] Filtering " << stat << "::" << param_name << "::" << filtername << " " << tmp[kk].date.toString(Date::ISO_TZ) << " [" << tmp[kk].date.toString(Date::ISO_WEEK) << "]\n";
						}
					}
					#endif
					if ((jj+1) != nr_of_filters){//after the last filter not necessary
						for (size_t kk=0; kk<ovec[ii].size(); kk++){
							tmp[kk](param) = ovec[ii][kk](param);
						}
					}
				} else {
					ostringstream ss;
					ss << "The filter \"" << (*filter_stack[jj]).getName() << "\" received " << tmp.size();
					ss << " timestamps and returned " << ovec[ii].size() << " timestamps!";
					throw IndexOutOfBoundsException(ss.str(), AT);
				}
			}

			if (!appliedFilter) //if not a single filter was applied
				ovec[ii] = ivec[ii]; //just copy input to output
		} else {
			ovec[ii] = ivec[ii]; //just copy input to output
		}
	}
}

const std::string ProcessingStack::toString() const
{
	std::ostringstream os;
	//os << "<ProcessingStack>";
	os << setw(10) << param_name << "::";

	for (size_t ii=0; ii<filter_stack.size(); ii++) {
		os << setw(10) << (*filter_stack[ii]).toString();
	}

	//os << "</ProcessingStack>";
	os << "\n";
	return os.str();
}

} //end namespace
