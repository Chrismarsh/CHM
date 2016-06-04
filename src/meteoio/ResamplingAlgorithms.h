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
#ifndef __RESAMPLINGALGORITHMS_H__
#define __RESAMPLINGALGORITHMS_H__

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/StationData.h>
#include <meteoio/meteoStats/libinterpol1D.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>

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
 * Two algorithms for the resampling are implemented:
 * - none: do not perform resampling, see NoResampling
 * - nearest:  nearest neighbor data resampling, see NearestNeighbour
 * - linear: linear data resampling, see LinearResampling
 * - accumulate: data re-accumulation as suitable for precipitations, see Accumulate
 * - daily_solar: generate solar radiation (ISWR or RSWR) from daily sums, see Daily_solar
 * - daily_avg: generate a sinusoidal variation around the measurement taken as daily average and of a given amplitude, see DailyAverage
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
};

class ResamplingAlgorithmsFactory {
	public:
		static ResamplingAlgorithms* getAlgorithm(const std::string& i_algoname, const std::string& parname, const double& window_size, const std::vector<std::string>& vecArgs);
};

/**********************************************************************************
 * The following functions are implementations of different resampling algorithms *
 **********************************************************************************/

/**
 * @brief No resampling: do not resample parameter but keep original sampling rate
 * It is enabled either with the "none" or "no" key.
 * @code
 * [Interpolations1D]
 * TA::resample = none
 * @endcode
 */
class NoResampling : public ResamplingAlgorithms {
	public:
		NoResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs);

		void resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
		std::string toString() const;
};

/**
 * @brief Nearest Neighbour data resampling
 * Find the nearest neighbour of a desired data point that is not IOUtils::nodata and copy that value into the desired data point
 *        - If the data point itself is not IOUtils::nodata, nothing needs to be done
 *        - If two points have the same distance from the data point to be resampled, calculate mean and return it
 *        - if the argument extrapolate is provided, points within WINDOW_SIZE seconds of only one valid point will receive the value of this point
 * The window size can be specified as argument but must appear in first position.
 * @code
 * [Interpolations1D]
 * TA::resample   = nearest
 * TA::nearest = 86400 extrapolate
 * @endcode
 */
class NearestNeighbour : public ResamplingAlgorithms {
	public:
		NearestNeighbour(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs);

		void resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
		std::string toString() const;
	private:
		bool extrapolate;
};

/**
 * @brief Linear data resampling: If a point is requested that is in between two input data points,
 *        the requested value is automatically calculated using a linear interpolation. Furthermore
 *        if the argument extrapolate is provided there will be an attempt made to extrapolate the
 *        point if the interpolation fails, by solving the line equation y = kx + d
 * The window size can be specified as argument but must appear in first position.
 * @code
 * [Interpolations1D]
 * TA::resample = linear
 * TA::linear   = extrapolate
 * @endcode
 */
class LinearResampling : public ResamplingAlgorithms {
	public:
		LinearResampling(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs);

		void resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
		std::string toString() const;
	private:
		bool extrapolate;
};

/**
 * @brief Accumulation over a user given period.
 * The input data is accumulated over a given time interval (given as filter argument, in seconds).
 * This is for example needed for converting rain gauges measurements read every 10 minutes to
 * hourly precipitation measurements. Remarks:
 * - the accumulation period has to be provided as an argument (in seconds)
 * - if giving the argument "strict", nodatas will propagate (ie. a single nodata in the input will force the re-accumulated value to be nodata). By default, all valid values are aggregated and only pure nodata intervals produce a nodata in the output.
 * @code
 * PSUM::resample   = accumulate
 * PSUM::accumulate = 3600
 * @endcode
 */
class Accumulate : public ResamplingAlgorithms {
	public:
		Accumulate(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs);

		void resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
		std::string toString() const;
	private:
		static size_t findStartOfPeriod(const std::vector<MeteoData>& vecM, const size_t& index, const Date& dateStart);
		double easySampling(const std::vector<MeteoData>& vecM, const size_t& paramindex, const size_t& /*index*/, const size_t& start_idx, const Date& dateStart, const Date& resampling_date) const;
		double complexSampling(const std::vector<MeteoData>& vecM, const size_t& paramindex, const size_t& index, const size_t& start_idx, const Date& dateStart, const Date& resampling_date) const;

		double accumulate_period; //internally, in julian days
		bool strict;
};

/**
 * @brief Generate solar radiation out of daily sums.
 * Daily sums of solar radiation (once, per day, any time during the day). Data provided at midnight is considered to belong to the day that just finished)
 * are compared to the potential radiation, leading to an atmospheric loss factor.
 * This loss factor is then applied to the potential solar radiation calculated at the requested time.
 * When using this algorithm for RSWR, an albedo is required. A default value of 0.5 is used. If the snow height is available and greater than a 10cm threshold,
 * a snow albedo is used. Below this threshold, a soil albedo is used.
 * @code
 * ISWR::resample   = daily_solar
 * @endcode
 */
class Daily_solar : public ResamplingAlgorithms {
	public:
		Daily_solar(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs);

		void resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
		std::string toString() const;
	private:
		double getSolarInterpol(const Date& resampling_date, const size_t& stat_idx) const;
		double compRadiation(const double& lat, const double& lon, const double& alt, const double& HS, const size_t& stat_idx);
		size_t getStationIndex(const std::string& key);

		std::vector< std::vector<double> > radiation;
		std::vector<std::string> station_index;
		std::vector<Date> dateStart, dateEnd;
		std::vector<double> loss_factor;

		static const double soil_albedo, snow_albedo, snow_thresh;
		static const size_t samples_per_day;
};

/**
 * @brief Generate daily variations of a given amplitude around a single daily average.
 * The paremeter to be interpolated is assumed to be a daily average and a sinusoidal variation of the
 * amplitude given as argument will be generated (it is also possible to provide the "phase" or the
 * fraction of the day when the minimum is reached). If data bearing the same name followed by "_MIN" or "_MAX"
 * exist, there is no need to provide an amplitude as they will be used instead (but if the amplitude is provided, it
 * will be used as a fallback when no min or max is available).
 *
 * @code
 * [Interpolations1D]
 * TA::resample = daily_avg
 * TA::daily_avg = 5 .25                ;assume that TA varies +/- 5K around its average during the day and reaches its minimum at 6am
 * @endcode
 * @note If both the average (the parameter itself in the data set),
 * min and max are provided, an error message will be returned.
 */
class DailyAverage : public ResamplingAlgorithms {
	public:
		DailyAverage(const std::string& i_algoname, const std::string& i_parname, const double& dflt_window_size, const std::vector<std::string>& vecArgs);

		void resample(const size_t& index, const ResamplingPosition& position, const size_t& paramindex,
		              const std::vector<MeteoData>& vecM, MeteoData& md);
		std::string toString() const;
	private:
		double getValue(const std::vector<MeteoData>& vecM, const size_t& paramindex, const size_t& index, const Date& dayStart, const double& frac_day) const;
		double range, phase;
};

} //end namespace
#endif
