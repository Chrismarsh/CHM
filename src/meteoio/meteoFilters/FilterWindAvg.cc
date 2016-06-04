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
#include <cmath>
#include <meteoio/meteoFilters/FilterWindAvg.h>
#include <meteoio/meteoLaws/Meteoconst.h>

using namespace std;

namespace mio {

FilterWindAvg::FilterWindAvg(const std::vector<std::string>& vec_args, const std::string& name) : WindowedFilter(name)
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise:
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void FilterWindAvg::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                            std::vector<MeteoData>& ovec)
{
	if(param!=MeteoData::VW && param!=MeteoData::DW) {
		ostringstream ss;
		ss << "Can not use " << getName() << " processing on " << MeteoData::getParameterName(param);
		throw InvalidArgumentException(ss.str(), AT);
	}

	ovec = ivec;

	for (size_t ii=0; ii<ovec.size(); ii++){ //for every element in ivec, get a window
		double& value = ovec[ii](param);

		size_t start, end;
		if( get_window_specs(ii, ivec, start, end) ) {
			value = calc_avg(ivec, param, start, end);
		} if(!is_soft) value = IOUtils::nodata;
	}
}

double FilterWindAvg::calc_avg(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end)
{
	//calculate ve and vn
	double ve=0.0, vn=0.0;
	size_t count=0;
	for (size_t ii=start; ii<=end; ii++) {
		const double VW = ivec[ii](MeteoData::VW);
		const double DW = ivec[ii](MeteoData::DW);
		if(VW!=IOUtils::nodata && DW!=IOUtils::nodata) {
			ve += VW * sin(DW*Cst::to_rad);
			vn += VW * cos(DW*Cst::to_rad);
			count++;
		}
	}

	if(count==0) return IOUtils::nodata;

	ve /= static_cast<double>(count);
	vn /= static_cast<double>(count);

	if(param==MeteoData::VW) {
		const double meanspeed = sqrt(ve*ve + vn*vn);
		return meanspeed;
	} else {
		const double meandirection = fmod( atan2(ve,vn) * Cst::to_deg + 360. , 360.); // turn into degrees [0;360)
		return meandirection;
	}
}

void FilterWindAvg::parse_args(std::vector<std::string> vec_args)
{
	vector<double> filter_args;

	if (vec_args.size() > 2){
		is_soft = ProcessingBlock::is_soft(vec_args);
	}

	if (vec_args.size() > 2)
		centering = (WindowedFilter::Centering)WindowedFilter::get_centering(vec_args);

	convert_args(2, 2, vec_args, filter_args);

	if ((filter_args[0] < 1) || (filter_args[1] < 0)){
		throw InvalidArgumentException("Invalid window size configuration for filter " + getName(), AT);
	}

	min_data_points = (unsigned int)floor(filter_args[0]);
	min_time_span = Duration(filter_args[1] / 86400.0, 0.);
}

}
