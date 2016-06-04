/***********************************************************************************/
/*  Copyright 2013 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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

#ifndef __DATAGENERATOR_H__
#define __DATAGENERATOR_H__

#include <meteoio/Config.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/GeneratorAlgorithms.h>

#include <vector>
#include <map>

namespace mio {

/**
 * @page dev_DataGenerator How to write a data generator
 * Once the data has been read, filtered and resampled, it can be that some data points are still missing.
 * These are either a few isolated periods (a sensor was not functioning) that are too large for performing
 * a statistical temporal interpolation or that a meteorological parameter was not even measured. In such a case,
 * we generate data, generally relying on some parametrization using other meteorological parameters. In a few
 * cases, even fully arbitrary data might be helpful (replacing missing value by a given constant so a model can
 * run over the data gap).
 *
 * @section structure_DataGenerator Structure
 * The selection of which data generator to use at any given time step, for a given parameter is
 * performed by the DataGenerator class. This class acts as an interface, presenting a higher level view to the
 * caller. The data generators themselves derive from the GeneratorAlgorithm class that standardizes their
 * public API. An object factory creates the generator during intialization (keeping all constructed generators
 * in a vector during the whole life time of the DataGenerator object), based on the strings contained in the user's
 * io.ini configuration file.
 *
 * The API also defines two public "generate" methods, taking a meteorological parameter index (see MeteoData) and either
 * a set of meteo data for one station and at one point in time or a meteo time series for one station.
 * These methods walk through the meteo data looking for nodata values for the requested meteo parameter index.
 * If the generator could successfully generate data for <b>all</b> the nodata values it found, it returns <i>true</i>,
 * <i>false</i> otherwise. If <i>false</i> was returned, the DataGenerator object that manages the process would
 * call the next data generator, <b>in the order that was declared by the user</b>. For a given meteo parameter, the
 * whole process stops as soon as a <i>true</i> is returned or there are no more data generators to try
 * (as declared by the user in his configuration file).
 *
 * @section implementation_DataGenerator Implementation
 * It is therefore necessary to create in GeneratorAlgorithms.cc (and declared in the .h) a new class,
 * nammed after the generator that will be implemented and inheriting GeneratorAlgorithm. Three methods need
 * to be implemented:
 * - the constructor with (const std::vector<std::string>& vecArgs, const std::string& i_algo)
 * - bool generate(const size_t& param, MeteoData& md)
 * - bool generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
 *
 * The constructor is responsible for parsing the arguments as a vector of strings and saving its own name internally, for
 * error messages, warnings, etc. It should set all internal variables it sees fit according to the parsed arguments. The
 * goal is <i>to not</i> do any parsing anywhere else (for performances reasons).
 *
 * The <i>generate(const size_t& param, MeteoData& md)</i> method compares <i>md(param)</i> with <i>IOUtils::nodata</i> and replaces
 * it by its generated value if necessary. It returns <i>true</i> if no further processing is needed
 * (ie. no replacement was needed or the replacement could be done) or <i>false</i> otherwise.
 *
 * The <i>generate(const size_t& param, std::vector<MeteoData>& vecMeteo)</i> method compares
 * <i>vecMeteo[ii](param)</i> with <i>IOUtils::nodata</i> for each timestamp in the vector and tries to generate data when necessary.
 * If all missing data points could be generated (or if no data point required to be generated), it returns <i>true</i>,
 * and <i>false</i> otherwise.
 *
 * Finally, a new entry must be added in the object factory GeneratorAlgorithmFactory::getAlgorithm method at the top of file
 * GeneratorAlgorithms.cc.
 *
 * @section doc_DataGenerator Documentation
 * The newly added data generator must be added to the list of available algorithms in
 * GeneratorAlgorithms.h with a proper description. Its class must be properly documented, similarly to the other data
 * generators. An example can also be given in the example section of the same file.
 * Please feel free to add necessary bibliographic references to the bibliographic section!
 */

/**
 * @class DataGenerator
 * @brief A class to generate meteo data from user-selected models or parametrizations.
 * This class sits in between the actual implementation of the various methods and the IOManager in
 * order to offer some high level interface. It basically reads the arguments and creates the objects for
 * the various data generators in its constructor and loop through the parameters and stations when called to fill the data.
 *
 * @ingroup meteoLaws
 * @author Mathias Bavay
 * @date   2013-03-20
 */

class DataGenerator {
	public:
		DataGenerator(const Config& cfg);
		DataGenerator(const DataGenerator& c) : mapGenerators(c.mapGenerators), mapCreators(c.mapCreators),
		                                        generators_defined(c.generators_defined), creators_defined(c.creators_defined) {}
		~DataGenerator();

		void fillMissing(METEO_SET& vecMeteo) const;
		void fillMissing(std::vector<METEO_SET>& vecVecMeteo) const;

		void createParameters(METEO_SET& vecMeteo) const;
		void createParameters(std::vector<METEO_SET>& vecVecMeteo) const;

		DataGenerator& operator=(const DataGenerator& source);

		const std::string toString() const;

	private:
		static void getParameters(const Config& cfg, const std::string& key_pattern, std::set<std::string>& set_parameters);
		static size_t getAlgorithmsForParameter(const Config& cfg, const std::string& key_pattern, const std::string& parname, std::vector<std::string>& vecAlgorithms);
		static size_t getArgumentsForAlgorithm(const Config& cfg, const std::string& parname,
		                                const std::string& algorithm,
		                                std::vector<std::string>& vecArgs);
		static void setAlgorithms(const Config& cfg, const std::string& key_pattern, std::map< std::string, std::vector<GeneratorAlgorithm*> > &mapAlgorithms);

		std::map< std::string, std::vector<GeneratorAlgorithm*> > mapGenerators; //per parameter data generators algorithms
		std::map< std::string, std::vector<GeneratorAlgorithm*> > mapCreators; //per parameter data creators algorithms
		bool generators_defined, creators_defined; //if true, there are some generators to run. if false, nothing to do
};

} //end namespace

#endif
