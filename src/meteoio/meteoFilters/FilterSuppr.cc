/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <ctime>
#include <cstdlib>

#include <meteoio/meteoFilters/FilterSuppr.h>

using namespace std;

namespace mio {

FilterSuppr::FilterSuppr(const std::vector<std::string>& vec_args, const std::string& name)
          : FilterBlock(name), range(IOUtils::nodata)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void FilterSuppr::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	
	if (range==IOUtils::nodata) { //remove all
		for (size_t ii=0; ii<ovec.size(); ii++){
			ovec[ii](param) = IOUtils::nodata;
		}
	} else { //only remove a given fraction
		const size_t set_size = ovec.size();
		const size_t nrRemove = static_cast<size_t>( round( (double)set_size*range ) );

		srand( static_cast<unsigned int>(time(NULL)) );
		size_t ii=1;
		while(ii<nrRemove) {
			const size_t idx = rand() % set_size;
			if (ivec[idx](param)!=IOUtils::nodata && ovec[idx](param)==IOUtils::nodata) continue; //the point was already removed
			
			ovec[idx](param)=IOUtils::nodata; //ie nodata points remain and are counted
			ii++;
		}
	}
}

void FilterSuppr::parse_args(std::vector<std::string> vec_args) 
{
	const size_t nrArgs = vec_args.size();

	if (nrArgs>1)
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);
	
	if (nrArgs==1) {
		if (!IOUtils::convertString(range, vec_args[0]))
			throw InvalidArgumentException("Invalid range \""+vec_args[0]+"\" specified for the "+getName()+" filter.", AT);
		if (range<0. || range>1.)
			throw InvalidArgumentException("Wrong range for filter " + getName() + ", it should be between 0 and 1", AT);
	}
}

} //end namespace
