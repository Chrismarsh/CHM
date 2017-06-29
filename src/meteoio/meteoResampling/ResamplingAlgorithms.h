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
#ifndef RESAMPLINGALGORITHMS_H
#define RESAMPLINGALGORITHMS_H

#include <meteoio/dataClasses/MeteoData.h>

#include <string>
#include <vector>

#ifdef _MSC_VER
	#pragma warning(disable:4512) //we don't need any = operator!
#endif

namespace mio {
 /**
 * @page resampling Resampling overview
 * The resampling infrastructure is described in ResamplingAlgorithms (for its API).
 * The goal of this page is to give an overview of the available resampling algorithms and their usage.
 *
 * @section resampling_section Resampling section
 * The resampling is specified for each parameter in the [Interpol1D] section. This section contains
 * a list of the various meteo parameters with their associated choice of resampling algorithm and
 * optional parameters. If a meteo parameter is not listed in this section, a linear resampling would be
 * assumed. An example of such section is given below:
 * @code
 * [Interpolations1D]
 * WINDOW_SIZE     = 86400
 * TA::resample    = linear
 *
 * RH::resample    = linear
 * RH::linear      = 172800
 *
 * VW::resample    = n_neighbor
 * VW::n_neighbor  = extrapolate
 *
 * PSUM::resample   = accumulate
 * PSUM::accumulate = 3600
 * @endcode
 *
 * Most of the resampling algorithms allow you to define per-meteo parameter and per-algorithm the WINDOW_SIZE. Otherwise, the section's WINDOW_SIZE is
 * used as default window size. This represents the biggest gap that can be interpolated (in seconds). Therefore if two valid points are less than
 * WINDOW_SIZE seconds apart, points in between will be interpolated. If they are further apart, all points in between will remain IOUtils::nodata.
 * If using the "extrapolate" optional argument, points at WINDOW_SIZE distance of only one valid point will be extrapolated, otherwise they will remain
 * IOUtils::nodata. Please keep in mind that allowing extrapolated values can lead to grossly out of range data: using the slope
 * between two hourly measurements to extrapolate a point 10 days ahead is obviously risky!
 *
 * By default, WINDOW_SIZE is set to 2 days. This key has a <b>potentially large impact on run time/performance</b>.
 *
 * @section algorithms_available Available Resampling Algorithms
 * Several algorithms for the resampling are implemented:
 * - none: do not perform resampling, see NoResampling
 * - nearest:  nearest neighbor data resampling, see NearestNeighbour
 * - linear: linear data resampling, see LinearResampling
 * - accumulate: data re-accumulation as suitable for precipitations, see Accumulate
 * - solar: resample solar radiation by interpolating an atmospheric loss factor, see Solar
 * - daily_solar: generate solar radiation (ISWR or RSWR) from daily sums, see Daily_solar
 * - daily_avg: generate a sinusoidal variation around the measurement taken as daily average and of a given amplitude, see DailyAverage
 * 
 * By default a linear resampling will be performed. It is possible to turn off all resampling by setting the *Enable_Resampling* key 
 * to *false* in the [Interpolations1D] section.
 */

/**
 * @class ResamplingAlgorithms
 * @brief Interface class for the temporal resampling algorithms
 * These models generate data points that are missing based on neighbouring points in a time series.
 *
 * @ingroup stats
 * @author Mathias Bavay - Thomas Egger
 * @date   2013-05-24
 */
class ResamplingAlgorithms {

	public:
		enum ResamplingPosition {
			exact_match,
			before,
			after,
			begin,
			end
		};

		ResamplingAlgorithms(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& /*vecArgs*/)
		                    : algo(i_algoname), parname(i_parname), window_size(dflt_window_size) {}

		virtual ~ResamplingAlgorithms() {}

		const std::string getAlgo() const {return algo;}

		virtual void resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md) = 0;

		virtual std::string toString() const = 0;

 	protected:
		static double partialAccumulateAtLeft(const std::vector<MeteoData>& vecM, const size_t& paramindex,
		                                      const size_t& pos, const Date& curr_date);
		static double partialAccumulateAtRight(const std::vector<MeteoData>& vecM, const size_t& paramindex,
		                                       const size_t& pos, const Date& curr_date);
		static void getNearestValidPts(const size_t& pos, const size_t& paramindex, const std::vector<MeteoData>& vecM, const Date& resampling_date,
		                               const double& window_size, size_t& indexP1, size_t& indexP2);
		static double linearInterpolation(const double& x1, const double& y1,
		                                  const double& x2, const double& y2, const double& x3);
		static Date getDailyStart(const Date& resampling_date);
		static size_t getDailyValue(const std::vector<MeteoData>& vecM, const size_t& paramindex, size_t pos, const Date& intervalStart, const Date& intervalEnd);

		const std::string algo, parname;
		double window_size;
		static const double soil_albedo, snow_albedo, snow_thresh; ///< These thresholds are used to handle solar radiation
};

class ResamplingAlgorithmsFactory {
	public:
		static ResamplingAlgorithms* getAlgorithm(const std::string& i_algoname, const std::string& parname, const double& window_size, const std::vector<std::string>& vecArgs);
};

} //end namespace
#endif
