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
#include <meteoio/meteoFilters/FilterMAD.h>
#include <meteoio/meteoStats/libinterpol1D.h>
#include <cmath>

using namespace std;

namespace mio {

FilterMAD::FilterMAD(const std::vector<std::string>& vec_args, const std::string& name) : WindowedFilter(name)
{
	parse_args(vec_args);

	//This is safe, but maybe too imprecise
	properties.time_before = min_time_span;
	properties.time_after  = min_time_span;
	properties.points_before = min_data_points;
	properties.points_after = min_data_points;
}

void FilterMAD::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	for (size_t ii=0; ii<ovec.size(); ii++){ //for every element in ivec, get a window
		double& value = ovec[ii](param);
		if (value==IOUtils::nodata) continue;

		size_t start, end;
		if ( get_window_specs(ii, ivec, start, end) ) {
			MAD_filter_point(ivec, param, start, end, value);
		} else if (!is_soft) value = IOUtils::nodata;
	}
}

void FilterMAD::MAD_filter_point(const std::vector<MeteoData>& ivec, const unsigned int& param, const size_t& start, const size_t& end, double &value)
{
	const double K = 1. / 0.6745;

	std::vector<double> data( end-start+1 );
	for (size_t ii=start; ii<=end; ii++) data[ii-start] = ivec[ii](param);

	//Calculate MAD
	const double median = Interpol1D::getMedian(data);
	const double mad    = Interpol1D::getMedianAverageDeviation(data);

	if ( median==IOUtils::nodata || mad==IOUtils::nodata ) return;

	const double sigma = mad * K;
	const double upper_lim = median + 3.*sigma;
	const double lower_lim = median - 3.*sigma;

	if ( (value>upper_lim) || (value<lower_lim) ) {
		value = IOUtils::nodata;
	}
}

void FilterMAD::parse_args(std::vector<std::string> vec_args)
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
