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

#include <meteoio/DataCreator.h>

using namespace std;

namespace mio {

DataCreator::DataCreator(const Config& cfg)
              : mapCreators()
{
	const std::string key_pattern( "::create" );
	std::set<std::string> set_of_used_parameters;
	getParameters(cfg, key_pattern, set_of_used_parameters);

	std::set<std::string>::const_iterator it;
	for (it = set_of_used_parameters.begin(); it != set_of_used_parameters.end(); ++it) {
		std::vector<std::string> tmpAlgorithms;
		const std::string parname( *it );
		const size_t nrOfAlgorithms = getAlgorithmsForParameter(cfg, key_pattern, parname, tmpAlgorithms);

		std::vector<GeneratorAlgorithm*> vecGenerators(nrOfAlgorithms);
		for (size_t jj=0; jj<nrOfAlgorithms; jj++) {
			std::vector<std::string> vecArgs;
			getArgumentsForAlgorithm(cfg, parname, tmpAlgorithms[jj], vecArgs);
			vecGenerators[jj] = GeneratorAlgorithmFactory::getAlgorithm( cfg, tmpAlgorithms[jj], vecArgs);
		}

		if (nrOfAlgorithms>0) {
			mapCreators[parname] = vecGenerators;
		}
	}
}

DataCreator::~DataCreator()
{ //we have to deallocate the memory allocated by "new GeneratorAlgorithm()"
	std::map< std::string, std::vector<GeneratorAlgorithm*> >::iterator it;

	for (it=mapCreators.begin(); it!=mapCreators.end(); ++it) {
		std::vector<GeneratorAlgorithm*> &vec( it->second );
		for (size_t ii=0; ii<vec.size(); ii++)
			delete vec[ii];
	}
}

DataCreator& DataCreator::operator=(const DataCreator& source)
{
	if (this != &source) {
		mapCreators = source.mapCreators;
	}
	return *this;
}

/**
 * @brief create new parameters from parametrizations
 * This relies on data creators defined by the user for each meteo parameters.
 * This loops over the defined generators and stops as soon as all points
 * have been successfully created.
 * @param vecVecMeteo vector containing a timeserie for each station
 */
void DataCreator::createParameters(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (mapCreators.empty()) return; //no creators defined by the end user

	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator it;
	for (it=mapCreators.begin(); it!=mapCreators.end(); ++it) {
		const std::vector<GeneratorAlgorithm*> vecGenerators( it->second );

		for (size_t station=0; station<vecVecMeteo.size(); ++station) { //process this parameter on all stations
			//create the new parameter
			for (size_t ii=0; ii<vecVecMeteo[station].size(); ++ii) {
				vecVecMeteo[station][ii].addParameter( it->first );
			}
			const size_t param = vecVecMeteo[station][0].getParameterIndex(it->first);

			//fill the new parameter
			size_t jj=0;
			while (jj<vecGenerators.size() && vecGenerators[jj]->create(param, vecVecMeteo[station]) != true) jj++;
		}
	}
}

void DataCreator::getParameters(const Config& cfg, const std::string& key_pattern, std::set<std::string>& set_parameters)
{
	set_parameters.clear();
	std::vector<std::string> vec_keys;
	cfg.findKeys(vec_keys, key_pattern, "Input", true); //search anywhere in key

	for (size_t ii=0; ii<vec_keys.size(); ii++) {
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			const string tmp( vec_keys[ii].substr(0,found) );
			set_parameters.insert( IOUtils::strToUpper(tmp) );
		}
	}
}

// This function retrieves the user defined generator algorithms for
// parameter 'parname' by querying the Config object
size_t DataCreator::getAlgorithmsForParameter(const Config& cfg, const std::string& key_pattern, const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	vecAlgorithms.clear();

	std::vector<std::string> vecKeys;
	cfg.findKeys(vecKeys, parname+key_pattern, "Input");

	if (vecKeys.size() > 1)
		throw IOException("Multiple definitions of " + parname + "::generators in config file", AT);;

	if (vecKeys.empty())
		return 0;

	cfg.getValue(vecKeys.at(0), "Input", vecAlgorithms, IOUtils::nothrow);

	return vecAlgorithms.size();
}

size_t DataCreator::getArgumentsForAlgorithm(const Config& cfg,
                                               const std::string& parname,
                                               const std::string& algorithm,
                                               std::vector<std::string>& vecArgs)
{
	vecArgs.clear();
	cfg.getValue(parname+"::"+algorithm, "Input", vecArgs, IOUtils::nothrow);

	return vecArgs.size();
}

const std::string DataCreator::toString() const {
	std::ostringstream os;
	os << "<DataCreator>\n";

	os << "Creators defined: " << std::boolalpha << !mapCreators.empty() << std::noboolalpha << "\n";
	if (!mapCreators.empty()) {
		os << "User list of creators:\n";
		for (std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator iter = mapCreators.begin(); iter != mapCreators.end(); ++iter) {
			os << setw(10) << iter->first << " :: ";
			for (size_t jj=0; jj<iter->second.size(); jj++) {
				os << iter->second[jj]->getAlgo() << " ";
			}
			os << "\n";
		}
	}

	os << "</DataCreator>\n";
	return os.str();
}

} //namespace
