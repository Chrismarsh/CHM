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
#include <meteoio/meteoFilters/FilterBlock.h>

namespace mio {

FilterBlock::FilterBlock(const std::string& filter_name) : ProcessingBlock(filter_name) {

}

FilterBlock::~FilterBlock() {}

void FilterBlock::extract_dbl_vector(const unsigned int& param, const std::vector<MeteoData>& ivec,
                                     std::vector<double>& ovec)
{
	ovec.resize( ivec.size() );
	for(size_t ii=0; ii<ivec.size(); ii++) {
		ovec[ii] = ivec[ii](param);
	}
}

void FilterBlock::extract_dbl_vector(const unsigned int& param, const std::vector<const MeteoData*>& ivec,
                                     std::vector<double>& ovec)
{
	ovec.resize( ivec.size() );
	for(size_t ii=0; ii<ivec.size(); ii++) {
		ovec[ii] =  (*ivec[ii])(param);
	}
}

} //end namespace
