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

#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H

#include <meteoio/Config.h>
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataGenerators/GeneratorAlgorithms.h>

#include <vector>
#include <map>
#include <set>

namespace mio {

/**
 * @class DataGenerator
 * @brief A class to generate meteo data from user-selected models or parametrizations.
 * This class sits in between the actual implementation of the various methods and the IOManager in
 * order to offer some high level interface. It basically reads the arguments and creates the objects for
 * the various data generators in its constructor and loop through the parameters and stations when called to fill the data.
 *
 * @ingroup meteoLaws
 * @author Mathias Bavay
 */

class DataGenerator {
	public:
		DataGenerator(const Config& cfg);
		DataGenerator(const DataGenerator& c) : mapGenerators(c.mapGenerators) {}
		~DataGenerator();

		void fillMissing(METEO_SET& vecMeteo) const;
		void fillMissing(std::vector<METEO_SET>& vecVecMeteo) const;

		DataGenerator& operator=(const DataGenerator& source);

		const std::string toString() const;

	private:
		static void getParameters(const Config& cfg, const std::string& key_pattern, std::set<std::string>& set_parameters);
		static size_t getAlgorithmsForParameter(const Config& cfg, const std::string& key_pattern, const std::string& parname, std::vector<std::string>& vecAlgorithms);
		static size_t getArgumentsForAlgorithm(const Config& cfg, const std::string& parname,
		                                const std::string& algorithm,
		                                std::vector<std::string>& vecArgs);

		std::map< std::string, std::vector<GeneratorAlgorithm*> > mapGenerators; //per parameter data generators algorithms
};

} //end namespace

#endif
