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

#include <meteoio/DataGenerator.h>

using namespace std;

namespace mio {

DataGenerator::DataGenerator(const Config& cfg)
              : mapGenerators()
{
	const std::string key_pattern( "::generators" );
	std::set<std::string> set_of_used_parameters;
	getParameters(cfg, key_pattern, set_of_used_parameters);

	std::set<string>::const_iterator it;
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
			mapGenerators[parname] = vecGenerators;
		}
	}
}

DataGenerator::~DataGenerator()
{ //we have to deallocate the memory allocated by "new GeneratorAlgorithm()"
	std::map< std::string, std::vector<GeneratorAlgorithm*> >::iterator it;

	for (it=mapGenerators.begin(); it!=mapGenerators.end(); ++it) {
		std::vector<GeneratorAlgorithm*> &vec( it->second );
		for (size_t ii=0; ii<vec.size(); ii++)
			delete vec[ii];
	}
}

DataGenerator& DataGenerator::operator=(const DataGenerator& source)
{
	if (this != &source) {
		mapGenerators = source.mapGenerators;
	}
	return *this;
}

/**
 * @brief generate data to fill missing data points.
 * This relies on data generators defined by the user for each meteo parameters.
 * This loops over the defined generators and stops as soon as all missing points
 * have been successfully replaced.
 * @param vecMeteo vector containing one point for each station
 */
void DataGenerator::fillMissing(METEO_SET& vecMeteo) const
{
	if (mapGenerators.empty()) return; //no generators defined by the end user

	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator it;
	for (it=mapGenerators.begin(); it!=mapGenerators.end(); ++it) {
		const std::vector<GeneratorAlgorithm*> vecGenerators( it->second );

		for (size_t station=0; station<vecMeteo.size(); ++station) { //process this parameter on all stations
			const size_t param = vecMeteo[station].getParameterIndex(it->first);
			if (param==IOUtils::npos) continue;

			#ifdef DATA_QA
			const double old_val = vecMeteo[station](param);
			const std::string statName( vecMeteo[station].meta.getStationName() );
			const std::string statID( vecMeteo[station].meta.getStationID() );
			const std::string stat = (!statID.empty())? statID : statName;
			#endif

			bool status = false;
			size_t jj=0;
			while (jj<vecGenerators.size() && status != true) { //loop over the generators
				status = vecGenerators[jj]->generate(param, vecMeteo[station]);
				jj++;
				#ifdef DATA_QA
				if (vecMeteo[station](param) != old_val) {
					const std::string parname( it->first );
					const std::string algo_name( vecGenerators[jj-1]->getAlgo() );
					const Date date( vecMeteo[station].date );
					cout << "[DATA_QA] Generating " << stat << "::" << parname << "::" << algo_name << " " << date.toString(Date::ISO_TZ) << " [" << date.toString(Date::ISO_WEEK) << "]\n";
				}
				#endif
			}
		}
	}
}

/**
 * @brief generate data to fill missing data points.
 * This relies on data generators defined by the user for each meteo parameters.
 * This loops over the defined generators and stops as soon as all missing points
 * have been successfully replaced.
 * @param vecVecMeteo vector containing a timeserie for each station
 */
void DataGenerator::fillMissing(std::vector<METEO_SET>& vecVecMeteo) const
{
	if (mapGenerators.empty()) return; //no generators defined by the end user

	std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator it;
	for (it=mapGenerators.begin(); it!=mapGenerators.end(); ++it) {
		const std::vector<GeneratorAlgorithm*> vecGenerators( it->second );

		for (size_t station=0; station<vecVecMeteo.size(); ++station) { //process this parameter on all stations
			const size_t param = vecVecMeteo[station][0].getParameterIndex(it->first);
			if (param==IOUtils::npos) continue;

			#ifdef DATA_QA
			const METEO_SET old_val = vecVecMeteo[station];
			const std::string statName( old_val[0].meta.getStationName() );
			const std::string statID( old_val[0].meta.getStationID() );
			const std::string stat = (!statID.empty())? statID : statName;
			#endif

			bool status = false;
			size_t jj=0;
			while (jj<vecGenerators.size() && status != true) { //loop over the generators
				status = vecGenerators[jj]->create(param, vecVecMeteo[station]);
				jj++;
				#ifdef DATA_QA
				const std::string parname( it->first );
				const std::string algo_name( vecGenerators[jj-1]->getAlgo() );
				for (size_t kk=0; kk<old_val.size(); kk++) {
					if (old_val[kk](param) != vecVecMeteo[station][kk](param)) {
						cout << "[DATA_QA] Generating " << stat << "::" << parname << "::" << algo_name << " " << old_val[kk].date.toString(Date::ISO_TZ) << "\n";
					}
				}
				#endif
			}
		}
	}
}

void DataGenerator::getParameters(const Config& cfg, const std::string& key_pattern, std::set<std::string>& set_parameters)
{
	set_parameters.clear();
	std::vector<std::string> vec_keys;
	cfg.findKeys(vec_keys, key_pattern, "Generators", true); //search anywhere in key

	for (size_t ii=0; ii<vec_keys.size(); ii++) {
		const size_t found = vec_keys[ii].find_first_of(":");
		if (found != std::string::npos){
			const std::string tmp( vec_keys[ii].substr(0,found) );
			set_parameters.insert( IOUtils::strToUpper(tmp) );
		}
	}
}

// This function retrieves the user defined generator algorithms for
// parameter 'parname' by querying the Config object
size_t DataGenerator::getAlgorithmsForParameter(const Config& cfg, const std::string& key_pattern, const std::string& parname, std::vector<std::string>& vecAlgorithms)
{
	vecAlgorithms.clear();

	std::vector<std::string> vecKeys;
	cfg.findKeys(vecKeys, parname+key_pattern, "Generators");

	if (vecKeys.size() > 1) throw IOException("Multiple definitions of " + parname + "::generators in config file", AT);;
	if (vecKeys.empty()) return 0;

	cfg.getValue(vecKeys.at(0), "Generators", vecAlgorithms, IOUtils::nothrow);
	return vecAlgorithms.size();
}

size_t DataGenerator::getArgumentsForAlgorithm(const Config& cfg,
                                               const std::string& parname,
                                               const std::string& algorithm,
                                               std::vector<std::string>& vecArgs)
{
	vecArgs.clear();
	cfg.getValue(parname+"::"+algorithm, "Generators", vecArgs, IOUtils::nothrow);

	return vecArgs.size();
}

const std::string DataGenerator::toString() const {
	std::ostringstream os;
	os << "<DataGenerator>\n";
	os << "Generators defined: " << std::boolalpha << !mapGenerators.empty() << std::noboolalpha << "\n";
	if (!mapGenerators.empty()) {
		os << "User list of generators:\n";
		for (std::map< std::string, std::vector<GeneratorAlgorithm*> >::const_iterator iter = mapGenerators.begin(); iter != mapGenerators.end(); ++iter) {
			os << setw(10) << iter->first << " :: ";
			for (size_t jj=0; jj<iter->second.size(); jj++) {
				os << iter->second[jj]->getAlgo() << " ";
			}
			os << "\n";
		}
	}

	os << "</DataGenerator>\n";
	return os.str();
}

} //namespace
