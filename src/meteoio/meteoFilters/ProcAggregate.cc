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
#include <meteoio/meteoFilters/ProcAggregate.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <meteoio/meteoLaws/Meteoconst.h>

using namespace std;

namespace mio {

ProcAggregate::ProcAggregate(const std::vector<std::string>& vec_args, const std::string& name) 
              : WindowedFilter(name), type(mean_agg)
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise:
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void ProcAggregate::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                            std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	
	for (size_t ii=0; ii<ovec.size(); ii++){ //for every element in ivec, get a window
		double& value = ovec[ii](param);
		size_t start, end;
		if ( get_window_specs(ii, ivec, start, end) ) {
			switch (type) {
				case min_agg:
					value = calc_min(ivec, param, start, end); break;
				case max_agg:
					value = calc_max(ivec, param, start, end); break;
				case mean_agg:
					value = calc_mean(ivec, param, start, end); break;
				case median_agg:
					value = calc_median(ivec, param, start, end); break;
				case wind_avg_agg:
					value = calc_wind_avg(ivec, param, start, end); break;
				default:
					throw UnknownValueException("Unknown aggregation algorithm selected!", AT);
			}
		} else {
			if (!is_soft) value = IOUtils::nodata;
		}
	}
}

double ProcAggregate::calc_min(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const
{
	double min = Cst::dbl_max;
	for (size_t ii=start; ii<=end; ii++){
		const double& value = ivec[ii](param);
		if (value!=IOUtils::nodata && value<min) min = value;
	}
	
	if (min == Cst::dbl_max) return IOUtils::nodata;
	return min;
}

double ProcAggregate::calc_max(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const
{
	double max = Cst::dbl_min;
	for (size_t ii=start; ii<=end; ii++){
		const double& value = ivec[ii](param);
		if (value!=IOUtils::nodata && value>max) max = value;
	}

	if (max == Cst::dbl_min) return IOUtils::nodata;
	return max;
}

double ProcAggregate::calc_mean(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const
{
	double sum = 0;
	size_t counter = 0;
	for (size_t ii=start; ii<=end; ii++){
		const double& value = ivec[ii](param);
		if (value != IOUtils::nodata){
			sum += value;
			counter++;
		}
	}

	if (counter == 0) return IOUtils::nodata;
	return (sum / (double)counter);
}

double ProcAggregate::calc_median(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const
{
	vector<double> vecTemp;
	vecTemp.reserve( end-start+1 );
	for (size_t ii=start; ii<=end; ii++){ //get rid of nodata elements
		const double& value = ivec[ii](param);
		if (value != IOUtils::nodata)
			vecTemp.push_back(value);
	}

	return Interpol1D::getMedian(vecTemp, false);
}

double ProcAggregate::calc_wind_avg(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end) const
{
	//calculate ve and vn
	double ve=0.0, vn=0.0;
	size_t count=0;
	for (size_t ii=start; ii<=end; ii++) {
		const double VW = ivec[ii](MeteoData::VW);
		const double DW = ivec[ii](MeteoData::DW);
		if (VW!=IOUtils::nodata && DW!=IOUtils::nodata) {
			ve += VW * sin(DW*Cst::to_rad);
			vn += VW * cos(DW*Cst::to_rad);
			count++;
		}
	}

	if (count==0) return IOUtils::nodata;

	ve /= static_cast<double>(count);
	vn /= static_cast<double>(count);

	if (param==MeteoData::VW) {
		const double meanspeed = sqrt(ve*ve + vn*vn);
		return meanspeed;
	} else {
		const double meandirection = fmod( atan2(ve,vn) * Cst::to_deg + 360. , 360.); // turn into degrees [0;360)
		return meandirection;
	}
}

void ProcAggregate::parse_args(std::vector<std::string> vec_args)
{
	if (vec_args.empty())
		throw InvalidArgumentException("Invalid number of arguments for filter " + getName(), AT);
	
	const string str_type = IOUtils::strToUpper( vec_args[0] );
	if (str_type=="MIN") type=min_agg;
	else if (str_type=="MAX") type=max_agg;
	else if (str_type=="MEAN") type=mean_agg;
	else if (str_type=="MEDIAN") type=median_agg;
	else if (str_type=="WIND_AVG") type=wind_avg_agg;
	else
		throw InvalidArgumentException("Unknown type '"+str_type+"' for filter " + getName(), AT);
	vec_args.erase( vec_args.begin() );

	if (vec_args.size() > 2)
		is_soft = ProcessingBlock::is_soft(vec_args);
	if (vec_args.size() > 2)
		centering = (WindowedFilter::Centering)WindowedFilter::get_centering(vec_args);

	vector<double> filter_args;
	convert_args(2, 2, vec_args, filter_args);
	if ((filter_args[0] < 1) || (filter_args[1] < 0)){
		throw InvalidArgumentException("Invalid window size configuration for filter " + getName(), AT);
	}

	min_data_points = (unsigned int)floor(filter_args[0]);
	min_time_span = Duration(filter_args[1] / 86400.0, 0.);
}

} //namespace
