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
#include <meteoio/meteoFilters/FilterRate.h>
#include <cmath>

using namespace std;

namespace mio {

FilterRate::FilterRate(const std::vector<std::string>& vec_args, const std::string& name)
           : FilterBlock(name), min_rate_of_change(0.), max_rate_of_change(0.)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::both; //for the rest: default values
}

void FilterRate::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                           std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	size_t last_good = IOUtils::npos;

	//Find first point that is not IOUtils::nodata
	for (size_t ii=0; ii<ovec.size(); ii++){
		if (ovec[ii](param) != IOUtils::nodata){
			last_good = ii;
			break;
		}
	}

	if (last_good == IOUtils::npos) //can not find a good point to start
		return;

	for (size_t ii=(last_good+1); ii<ovec.size(); ii++) {
		double& curr_value       = ovec[ii](param);
		const double& prev_value = ovec[last_good](param);
		const double curr_time   = ovec[ii].date.getJulian();
		const double prev_time   = ovec[last_good].date.getJulian();

		if (curr_value == IOUtils::nodata)
			continue;

		const double local_rate = (curr_value-prev_value) / ((curr_time-prev_time+1e-12)*24.*3600.); //per seconds

		if( local_rate>max_rate_of_change || local_rate<min_rate_of_change ) {
			curr_value = IOUtils::nodata;
		} else {
			last_good = ii;
		}
	}
}

void FilterRate::parse_args(const std::vector<std::string>& vec_args) {
	vector<double> filter_args;
	convert_args(1, 2, vec_args, filter_args);

	const size_t nb_args = filter_args.size();
	if (nb_args == 2) {
		min_rate_of_change = filter_args[0];
		max_rate_of_change = filter_args[1];
	} else if(nb_args == 1) {
		min_rate_of_change = -filter_args[0];
		max_rate_of_change = filter_args[0];
	} else
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName() + " - Please provide 1 or 2 arguments!", AT);

}

}
