/***********************************************************************************/
/*  Copyright 2014 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#ifndef __METEO2DINTERPOLATOR_H__
#define __METEO2DINTERPOLATOR_H__

#include <meteoio/TimeSeriesManager.h>
#include <meteoio/Config.h>
#include <meteoio/dataClasses/Buffer.h>
#include <meteoio/dataClasses/Date.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/DEMObject.h>
#include <meteoio/InterpolationAlgorithms.h>

#include <memory>
#include <vector>
#include <map>

namespace mio {

class InterpolationAlgorithm;

/**
 * @page dev_2Dinterpol How to write a spatial interpolation algorithm
 * Point measurements can be spatially interpolated by MeteoIO, through the use of interpolation
 * algorithms. The user will then be able to choose for each meteorological parameter which
 * interpolations could be applicable and the system will choose (through a heuristic) which one
 * to apply at each time step (depending on the conditions of the moment, like the number of measurements).
 *
 * @section structure_2Dinterpol Structure
 * The selection of which interpolation algorithm to use at any given time step, for a given parameter is
 * performed by the Meteo2DInterpolator class. This class provides the interface for the spatial
 * interpolations. The interpolation algorithms themselves derive from the
 * InterpolationAlgorithm class that standardizes their public API (which would not be used by anything
 * but a Meteo2DInterpolator object). It contains a getQualityRating() method that must return a quality
 * index between 0 (algorithm not applicable) and 1 (perfect result if using this algorithm). This is currently
 * only based on extremely simple heuristics, using general knowledge about the applicability of the various
 * spatial interpolation methods depending on some obvious factors (number of measurement points, etc). The
 * Meteo2DInterpolator object will call this method from all the algorithms listed by the user (in his io.ini
 * configuration file) and keep the one that gets the highest score for interpolating the current parameter
 * at the current time step. The interpolation is then done calling the algorithm's calculate method.
 *
 * @section implementation_2Dinterpol Implementation
 * It is therefore necessary to create in InterpolationAlgorithms.cc (and declared in the .h) a new class,
 * nammed after the algorithm that will be implemented and inheriting InterpolationAlgorithm. Three methods need
 * to be implemented (the constructor being inherited from InterpolationAlgorithm and automatically called
 * by an object factory):
 * - void initialize(const MeteoData::Parameters& in_param)
 * - double getQualityRating()
 * - void calculate(Grid2DObject& grid)
 *
 * The initialize method takes the meteorological parameter that will be interpolated and set the param
 * private member to it. It then computes the private member nrOfMeasurments that contains the number of
 * stations that have this meteorological parameter available by calling getData(param, vecData, vecMeta), which
 * also fills the vectors vecData and vecMeta with the available data (as double) and metadata (as StationData).
 * Custom data preparation can obviously be done in this method.
 *
 * The calculate method must properly erase and reste the grid that it receives before filling it. If necessary,
 * (as is the case for precipitation, relative humidity and snow height, for example) the grid can be checked for min/max by
 * calling checkMinMax() at the end of Meteo2DInterpolator::interpolate.It can also add extra information about the
 * interpolation process (such as a regression coefficient or error estimate) to the InterpolationAlgorithm::info
 * stringstream (which will be made available to external programs, such as GUIs).
 *
 * The new class and its associated end user key must be used and its constructor called in AlgorithmFactory::getAlgorithm.
 * It is recommended that any generic statistical
 * spatial processing be implemented as a static class in libinterpol2D.cc so that it could be reused by other
 * algorithms (see for example Interpol2D::IDW and IDWCore). In any case, proper doxygen documentation
 * must be written alongside the implementation.
 *
 * @section doc_2Dinterpol Documentation
 * The newly added interpolation algorithm must be added to the list of available algorithms in
 * InterpolationAlgorithms.h with a proper description. An example can also be given in the example section
 * of the same file. Please feel free to add necessary bibliographic references to the bibliographic section!
 *
*/

/**
 * @class Meteo2DInterpolator
 * @brief A class to spatially interpolate meteo parameters. For more, see \ref interpol2d
 *
 * @ingroup stats
 * @author Mathias Bavay and Thomas Egger
 * @date   2010-01-14
 */

class Meteo2DInterpolator {
	public:
		/**
		* @brief Constructor.
		*/
		Meteo2DInterpolator(const Config& i_cfg, TimeSeriesManager& i_tsmanager, GridsManager& i_gridsmanager);

		~Meteo2DInterpolator();

		///Keywords for virtual stations strategy
		typedef enum VSTATIONS_POLICY {
			VSTATIONS, ///< extract virtual stations as specified in the ini file
			DOWNSCALING, ///< extract all grid points from a provided grid
			SMART_DOWNSCALING ///< extract all relevant grid points from a provided grid
		} vstations_policy;

		/**
		 * @brief A generic function that can interpolate for any given MeteoData member variable
		 *
		 * @param date date for which to interpolate
		 * @param dem Digital Elevation Model on which to perform the interpolation
		 * @param meteoparam Any MeteoData member variable as specified in the
		 * 				 enum MeteoData::Parameters (e.g. MeteoData::TA)
		 * @param result A Grid2DObject that will be filled with the interpolated data
		 */
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result);

		/**
		 * @brief A generic function that can interpolate for any given MeteoData member variable
		 *
		 * @param date date for which to interpolate
		 * @param dem Digital Elevation Model on which to perform the interpolation
		 * @param meteoparam Any MeteoData member variable as specified in the
		 * 				 enum MeteoData::Parameters (e.g. MeteoData::TA)
		 * @param result A Grid2DObject that will be filled with the interpolated data
		 * @param InfoString some information about the interpolation process (useful for GUIs)
		 */
		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
		                 Grid2DObject& result, std::string& InfoString);

		void interpolate(const Date& date, const DEMObject& dem, const MeteoData::Parameters& meteoparam,
                            const std::vector<Coords>& in_coords, std::vector<double>& result, std::string& info_string);

		/**
		 * @brief Retrieve the arguments vector for a given interpolation algorithm
		 * @param param the meteorological parameter that is concerned
		 * @param algorithm the desired algorithm
		 * @param vecArgs a vector of strings containing the arguments
		 */
		size_t getArgumentsForAlgorithm(const std::string& param,
		                                const std::string& algorithm,
		                                std::vector<std::string>& vecArgs) const;

		/**
		 * @brief Compute point measurements from grids following a given computing strategy
		 * @param strategy sampling/computing strategy
		 * @param i_date when to compute the virtual stations
		 * @param vecMeteo a vector of meteodata for the configured virtual stations
		 */
		size_t getVirtualMeteoData(const vstations_policy& strategy, const Date& i_date, METEO_SET& vecMeteo);

		const std::string toString() const;

	private:
		static Config stripVirtualConfig(const Config& cfg);
		static void checkMinMax(const double& minval, const double& maxval, Grid2DObject& gridobj);
		static void check_projections(const DEMObject& dem, const std::vector<MeteoData>& vec_meteo);
		static size_t get_parameters(const Config& cfg, std::set<std::string>& set_parameters);
		static size_t getAlgorithmsForParameter(const Config& cfg, const std::string& parname, std::vector<std::string>& vecAlgorithms);

		size_t getVirtualStationsData(const Date& i_date, METEO_SET& vecMeteo);
		void setAlgorithms();
		void initVirtualStations();

		const Config& cfg; ///< Reference to Config object, initialized during construction
		TimeSeriesManager& tsmanager; ///< Reference to TimeSeriesManager object, used for callbacks, initialized during construction
		GridsManager& gridsmanager; ///< Reference to GridsManager object, used for callbacks, initialized during construction
		GridBuffer grid_buffer;

		std::map< std::string, std::vector<InterpolationAlgorithm*> > mapAlgorithms; //per parameter interpolation algorithms

		std::vector<size_t> v_params; ///< Parameters for virtual stations
		std::vector<Coords> v_coords; ///< Coordinates for virtual stations
		std::vector<StationData> v_stations; ///< metadata for virtual stations
		std::map<Date, METEO_SET > virtual_point_cache;  ///< stores already resampled virtual data points

		bool algorithms_ready; ///< Have the algorithms objects been constructed?
		bool use_full_dem; ///< use full dem for point-wise spatial interpolations
		bool downscaling; ///< Are we downscaling meteo grids instead of interpolating stations' data?
		bool virtual_stations; ///< compute the meteo values at virtual stations
};

} //end namespace

#endif
