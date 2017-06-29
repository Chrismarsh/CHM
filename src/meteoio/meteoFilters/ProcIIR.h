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
#ifndef PROCIIR_H
#define PROCIIR_H

#include <meteoio/meteoFilters/FilterBlock.h>
#include <vector>
#include <string>

namespace mio {

/**
 * @class  ProcIIR
 * @ingroup processing
 * @brief Infinite Impulse Response (IIR) filter.
 * This filter can either be used as a low pass or high pass filter. It is based on a Critically Damped, 2 poles filter (considering that
 * it is better avoid overshooting even at the cost of a gentler falloff). It is possible to use it as a Low Pass (**LP**) or High Pass (**HP**)
 *
 * The cutoff <b>period</b> (defined as the frequency at a -3dB gain) is given in seconds as argument. The phase is removed by
 * bidirectional filtering, ie. running the filter twice, first backward and then forward (this also squares the amplitude response). But
 * it is possible to disable this bidirectional filtering by adding the "single_pass" argument.
 *
 * @code
 * HS::filter1	= IIR
 * HS::arg1	= LP 10800 ;3 hours
 * @endcode
 *
 * To know more: http://unicorn.us.com/trading/allpolefilters.html and http://www.dspguide.com/ch19/4.htm.
 */

class ProcIIR : public ProcessingBlock {
	public:
		ProcIIR(const std::vector<std::string>& vec_args, const std::string& name);

		virtual void process(const unsigned int& param, const std::vector<MeteoData>& ivec,
		                     std::vector<MeteoData>& ovec);

	private:
		typedef enum IIR_TYPE {
			BUTTERWORTH,
			CRITICALLY_DAMPED,
			BESSEL
		} IIR_Type;

		static void getFilterParameters(const IIR_Type& type, const bool& isLowPass, const double& n, double &g, double &p, double &c);
		static double filterPoint(const double& raw_val, const double A[3], const double B[3], std::vector<double> &X, std::vector<double> &Y);
		void computeCoefficients(const double& fs, const double& f0, double A[3], double B[3]) const;

		void parse_args(std::vector<std::string> vec_args);

		double cutoff;
		double g, p, c; ///< filter definition: number of passes, polynomial coefficients, 3dB cutoff correction
		bool bidirectional, low_pass;
};

} //end namespace

#endif
