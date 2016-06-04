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
#include <meteoio/meteoFilters/ProcButterworth.h>
#include <meteoio/MathOptim.h>
#include <cmath>

using namespace std;

namespace mio {

ProcButterworth::ProcButterworth(const std::vector<std::string>& vec_args, const std::string& name)
                  : ProcessingBlock(name), X(3, IOUtils::nodata), Y(3, IOUtils::nodata), cutoff(0.)
{
	parse_args(vec_args);
	properties.points_before = 2;
	properties.stage = ProcessingProperties::first;
}

void ProcButterworth::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if(ivec.size()<2 || cutoff==0.) return;

	const double days = ivec.back().date.getJulian() - ivec.front().date.getJulian();
	const size_t nr_data_pts = ivec.size();
	const double sampling_rate = static_cast<double>(nr_data_pts-1) / (days*24.*3600.); //in Hz

	double A[3], B[3];
	computeCoefficients(sampling_rate, 1./cutoff, A, B);

	for (size_t ii=0; ii<ovec.size(); ++ii){
		const double& raw_val = ivec[ii](param);

		//propagate in X and Y
		X[2] = X[1]; X[1] = X[0]; X[0] = raw_val;
		Y[2] = Y[1]; Y[1] = Y[0]; Y[0] = raw_val; //Y[0] will be overwritten but in case of nodata we still propagate a value
		if(X[2]==IOUtils::nodata || X[1]==IOUtils::nodata || X[0]==IOUtils::nodata) continue;
		if(Y[2]==IOUtils::nodata || Y[1]==IOUtils::nodata) continue;

		Y[0] = A[0]*X[0] + A[1]*X[1] + A[2]*X[2] - ( B[1]*Y[1] + B[2]*Y[2] );
		ovec[ii](param) = Y[0];
	}

}

void ProcButterworth::computeCoefficients(const double& samplerate, const double& f_cutoff, double A[3], double B[3]) const
{
	const double QcRaw  = (2. * Cst::PI * f_cutoff) / samplerate; // Find cutoff frequency in [0..PI]
	const double QcWarp = tan(QcRaw); // Warp cutoff frequency
	const double gain = 1. / (1.+Cst::Sqrt2/QcWarp + 2./Optim::pow2(QcWarp));

	B[0] = 1.;
	B[1] = (2. - 2. * 2./Optim::pow2(QcWarp)) * gain;
	B[2] = (1. - Cst::Sqrt2/QcWarp + 2./Optim::pow2(QcWarp)) * gain;

	A[0] = 1. * gain;
	A[1] = 2. * gain;
	A[2] = 1. * gain;
}

void ProcButterworth::parse_args(std::vector<std::string> vec_args)
{
	vector<double> filter_args;
	convert_args(1, 1, vec_args, filter_args);

	if (filter_args.size() != 1)
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);

	cutoff = filter_args[0];
}

}
