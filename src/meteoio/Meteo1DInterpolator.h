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
#ifndef __METEO1DINTERPOLATOR_H__
#define __METEO1DINTERPOLATOR_H__

#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/StationData.h>
#include <meteoio/Config.h>
#include <meteoio/ResamplingAlgorithms.h>
#include <meteoio/meteoFilters/ProcessingBlock.h>

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <utility>

namespace mio {

/**
 * @class Meteo1DInterpolator
 * @brief A class that can resample MeteoData objects
 *
 * @ingroup stats
 * @author Thomas Egger
 * @date   2010-06-24
 */

class Meteo1DInterpolator {
	public:

		/**
		* @brief The default constructor
		* Set up the interpolation algorithm for each parameter
		* Init tasklist: a vector that holds one std::string for each parameter,
		*                representing the interpolation algorithm that will be executed
		*                for the respective parameter
		*                e.g. tasklist for TA: linear
		* taskargs:      a vector that holds the respective arguments for the algorithms
		*                as a std::vector<std::string>, so there can be multiple arguments
		*
		* @param[in] in_cfg Config object that holds the MeteoFilter configuration in the [Filters] section
		*/
		Meteo1DInterpolator(const Config& in_cfg);

		~Meteo1DInterpolator();

		/**
		 * @brief A function that executes all the resampling algorithms that have been setup in the constructor
		 * @param[in] date The requested date for a MeteoData object (to be resampled if not present)
		 * @param[in] vecM A vector of MeteoData where the new object will be inserted if not present
		 * @param[in] md new MeteoData element, filled with the resampled values
		 * @return true if successfull, false if no resampling was possible (no element created)
		 */
		bool resampleData(const Date& date, const std::vector<MeteoData>& vecM, MeteoData& md);

		void getWindowSize(ProcessingProperties& o_properties) const;

		Meteo1DInterpolator& operator=(const Meteo1DInterpolator&); ///<Assignement operator
		const std::string toString() const;

 	private:
		std::string getInterpolationForParameter(const std::string& parname, std::vector<std::string>& vecArguments) const;

		const Config& cfg;
		double window_size; ///< In seconds
		std::map< std::string, ResamplingAlgorithms* > mapAlgorithms; //per parameter interpolation algorithms
};
} //end namespace

#endif
