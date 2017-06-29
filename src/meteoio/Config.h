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
#ifndef CONFIGREADER_H
#define CONFIGREADER_H

#include <meteoio/IOUtils.h>
#include <meteoio/IOExceptions.h>

#include <string>
#include <sstream>
#include <map>
#include <vector>

namespace mio {

/**
 * @class Config
 * @brief A class that reads a key/value file. These files (typically named *.ini) follow the INI file format standard (see http://en.wikipedia.org/wiki/INI_file) and have the following structure:
 * - they consist of 0 or more explicitly indicated sections, which start with a sectionname enclosed in square brackets
 *   e.g. [General] or [Filter]
 * - within each section there are 0 or more key value pairs defined: KEY = VALUE
 * - in this implementation each key is unique within its section
 * - lines that start with a semicolon ';' or a hash sign '#' are ignored (regarded as comments)
 * - empty lines are ignored
 * - if there is no section name given in a file, the default section called "GENERAL" is assumed
 * - a VALUE for a KEY can consist of multiple whitespace separated values (e.g. MYNUMBERS = 17.77 -18.55 8888 99.99)
 * 
 * @anchor config_import
 * It is possible to import another ini file, by specifying as many of the keys listed below as necessary.
 *   Please not that in order to prevent circular dependencies, it is not possible to import the same file several times.
 *      - IMPORT_BEFORE = {file and path to import}. This must take place before any non-import
 *        key or section header. This imports the specified file before processing the current file, allowing
 *        to overwrite the imported parameters in the current configuration file.
 *      - IMPORT_AFTER = {file and path to import}. This can occur anywhere and imports the specified file
 *        after processing the current file, allowing to overwrite the local parameters by the imported parameters.
 *
 * @author Thomas Egger & Mathias Bavay
 * @date   2008-11-30
 */

class ConfigProxy;

class Config {
	public:
		/**
		 * @brief Empty constructor. The user MUST later one fill the internal key/value map object
		 */
		Config();

		virtual ~Config() {}

		/**
		 * @brief Main constructor. The file is parsed and a key/value map object is internally created
		 * @param[in] filename_in string representing the absolute filename of the key/value file
		 */
		Config(const std::string& filename_in);

		/**
		 * @brief Write the Config object to a file
		 * @param filename The filename including the path, e.g. "/tmp/test.ini"
		 */
		void write(const std::string& filename) const;

		/**
		 * @brief Add the content of a file to the internal key/value map object
		 * @param[in] filename_in string representing the absolute filename of the key/value file
		 */
		void addFile(const std::string& filename_in);

		/**
		 * @brief Add the content of the given command line to the internal key/value map object
		 * @param[in] cmd_line string representing the command line to be parsed for key/value pairs or switches
		 */
		void addCmdLine(const std::string& cmd_line);

		/**
		 * @brief Add a specific key/value pair to the internal key/value map object.
		 *        key and section are case insensitive
		 * @param[in] key string representing the key to be added
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[in] value string representing the matching value to be added
		 */
		void addKey(const std::string& key, const std::string& section, const std::string& value);

		/**
		 * @brief Delete a specific key/value pair from the internal map object, key/section are case insensitive
		 * @param[in] key string representing the key to be added
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		*/
		void deleteKey(const std::string& key, const std::string& section);

		/**
		 * @brief Delete keys matching a specific pattern from the internal map object, key/section are case insensitive
		 * @param[in] keymatch A string representing the beginning of a key to search for
		 * @param[in] section A string defining which section to search through (default: GENERAL)
		 * @param[in] anywhere Match substring anywhere in the key string (default=false, ie at the begining only)
		 * @code
		 *  Config cfg("io.ini");
		 *  cfg.deleteKeys("STATION", "Input");
		 * @endcode
		*/
		void deleteKeys(const std::string& keymatch, const std::string& section, const bool& anywhere=false);

		/**
		 * @brief Returns the filename that the Config object was constructed with.
		 * @return The absolute filename of the key/value file.
		 */
		std::string getSourceName() const;

		/**
		 * @brief Returns the directory where the root configuration file is (needed to resolv relative paths).
		 * @return The absolute path to the root config file (resolved for symlinks, relative paths, etc).
		 */
		std::string getConfigRootDir() const;

		/**
		 * @brief Return if a given key exists in a given section
		 * @param[in] key string representing the key to be searched
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @return true if the key exists
		 */
		bool keyExists(const std::string& key, const std::string& section) const;

		/**
		 * @brief Print the content of the Config object (usefull for debugging)
		 * The Config is bound by "<Config>" and "</Config>" on separate lines
		 */
		const std::string toString() const;

		friend std::iostream& operator<<(std::iostream& os, const Config& cfg);
		friend std::iostream& operator>>(std::iostream& is, Config& cfg);

		template <typename T> std::vector<T> getValue(const std::string& key, const std::string& section,
		                                              const IOUtils::ThrowOptions& opt=IOUtils::dothrow) const
		{
			std::vector<T> tmp;
			getValue(key, section, tmp, opt);
			return tmp;
		}

		/**
		 * @brief Template function to retrieve a vector of values of class T for a certain key
		 * @code
		 * algorithms = lsm linregres idw_kriging\n
		 * MYNUMBERS = 17.87 19.89 -19.89 +7777.007
		 * @endcode
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[out] vecT a variable of class vector<T> into which the values for the corresponding key are saved
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key,
		                                    std::vector<T>& vecT,
		                                    const IOUtils::ThrowOptions& opt=IOUtils::dothrow) const
		{
			getValue(key, "GENERAL", vecT, opt);
		}

		/**
		 * @brief Template function to retrieve a vector of values of class T for a certain key
		 * @code
		 * algorithms = lsm linregres idw_kriging\n
		 * MYNUMBERS = 17.87 19.89 -19.89 +7777.007
		 * @endcode
		 * @param[in] key std::string representing a KEY in the key/value file
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[out] vecT a variable of class vector<T> into which the values for the corresponding key are saved
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key, const std::string& section,
		                                    std::vector<T>& vecT, const IOUtils::ThrowOptions& opt=IOUtils::dothrow) const
		{
			vecT.clear();
			const std::string new_key( IOUtils::strToUpper(key) );
			const std::string new_section( IOUtils::strToUpper(section) );

			try {
				IOUtils::getValueForKey<T>(properties, new_section + "::" + new_key, vecT, opt);
			} catch(const std::exception&){
				throw UnknownValueException("[E] Error in "+sourcename+": no value for key "+new_section+"::"+new_key, AT);
			}
		}

		/**
		 * @ brief A function that allows to retrieve a value for a key as return parameter (vectors of values too)
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 * @return A value of type T
		 */
		ConfigProxy get(const std::string& key, const IOUtils::ThrowOptions& opt=IOUtils::dothrow) const;

		/**
		 * @ brief A function that allows to retrieve a value for a key as return parameter (vectors of values too)
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 * @return A value of type T
		 *
		 * Example Usage:
		 * @code
		 * Config cfg("io.ini");
		 * vector<int> = cfg.get("DEPTHS", "INPUT", IOUtils::nothrow);
		 * string mystr = cfg.get("PATH", "OUTPUT");
		 * @endcode
		 */
		ConfigProxy get(const std::string& key, const std::string& section, const IOUtils::ThrowOptions& opt=IOUtils::dothrow) const;

		/**
		 * @brief Template function to retrieve a value of class T for a certain key
		 * @param[in] key std::string representing a KEY in the key/value file (default section "GENERAL" is assumed)
		 * @param[out] t a variable of class T into which the value for the corresponding key is saved (e.g. double, int, std::string)
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key, T& t, const IOUtils::ThrowOptions& opt=IOUtils::dothrow) const
		{
			getValue(key, "GENERAL", t, opt);
		}

		/**
		 * @brief Template function to retrieve a value of class T for a certain key
		 * @param[in] key std::string representing a KEY in the key/value file
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[out] t a variable of class T into which the value for the corresponding key is saved (e.g. double, int, std::string)
		 * @param[in] opt indicating whether an exception should be raised, when key is not present
		 */
		template <typename T> void getValue(const std::string& key, const std::string& section, T& t,
                                              const IOUtils::ThrowOptions& opt=IOUtils::dothrow) const
		{
			const std::string new_key( IOUtils::strToUpper(key) );
			const std::string new_section( IOUtils::strToUpper(section) );

			try {
				IOUtils::getValueForKey<T>(properties, new_section + "::" + new_key, t, opt);
			} catch(const std::exception&){
				throw UnknownValueException("[E] Error in "+sourcename+": no value for key "+new_section+"::"+new_key, AT);
			}
		}

		/**
		 * @brief Template function to retrieve a vector of values of class T for a certain key pattern
		 * @param[in] keystart std::string representing a pattern for the key in the key/value file
		 * @param[in] section std::string representing a section name; the key has to be part of this section
		 * @param[out] vecT a vector of class T into which the values for the corresponding keys are saved
		 */
		template <typename T> void getValues(const std::string& keystart, const std::string& section, std::vector<T>& vecT) const
		{
			vecT.clear();
			std::vector< std::string > vecKeys;
			const std::string new_section( IOUtils::strToUpper(section) );
			const size_t nr_keys = findKeys(vecKeys, keystart, new_section);

			for (size_t ii=0; ii<nr_keys; ++ii) {
				const std::string full_key = new_section + "::" + vecKeys[ii];
				T tmp;
				try {
					IOUtils::getValueForKey<T>(properties, full_key, tmp, IOUtils::dothrow);
				} catch(const std::exception&){
					throw UnknownValueException("[E] Error in "+sourcename+" reading key "+full_key, AT);
				}
				vecT.push_back( tmp );
			}
		}

		template <typename T> void getValues(const std::string& keystart, const std::string& section, std::vector<T>& vecT, std::vector<std::string>& vecKeys) const
		{
			vecT.clear();
			const std::string new_section( IOUtils::strToUpper(section) );
			const size_t nr_keys = findKeys(vecKeys, keystart, new_section);

			for (size_t ii=0; ii<nr_keys; ++ii) {
				const std::string full_key = new_section + "::" + vecKeys[ii];
				T tmp;
				try {
					IOUtils::getValueForKey<T>(properties, full_key, tmp, IOUtils::dothrow);
				} catch(const std::exception&){
					throw UnknownValueException("[E] Error in "+sourcename+" reading key "+full_key, AT);
				}
				vecT.push_back( tmp );
			}
		}

		/**
		 * @brief Function that searches for a given string within the keys of section (default: GENERAL)
		 *         it returns the number of matches (partial matches are considered) and writes all the keys
		 *         into a vector\<string\> that is handed to the function as reference
		 * @param[out] vecResult A vector that will hold all keys that partially match keystart
		 * @param[in] keymatch A string representing the beginning of a key to search for
		 * @param[in] section A string defining which section to search through (default: GENERAL)
		 * @param[in] anywhere Match substring anywhere in the key string (default=false, ie at the begining only)
		 * @code
		 *  vector<string> myVec;
		 *  size_t nrOfMatches = cfg.findKeys(myVec, "TA::", "Filters");
		 * @endcode
		 */
		size_t findKeys(std::vector<std::string>& vecResult,
		               const std::string& keymatch, std::string section, const bool& anywhere=false) const;

	private:
		void parseCmdLine(const std::string& cmd_line);
		void parseFile(const std::string& filename);
		void parseLine(const unsigned int& linenr, std::vector<std::string> &import_after, bool &accept_import_before, std::string &line, std::string &section);
		std::string extract_section(std::string key) const;
		std::string clean_import_path(const std::string& in_path) const;

		std::map<std::string, std::string> properties; //Save key value pairs
		std::vector<std::string> imported; //list of files already imported (to avoid circular references)
		std::string sourcename; //description of the data source for the key/value pair
		std::string configRootDir; //directory of the root config file
		static const char* defaultSection;
}; //end class definition Config

class ConfigProxy {
	public:
		const Config& proxycfg;
		const std::string& key;
		const std::string& section;
		const IOUtils::ThrowOptions& opt;

		ConfigProxy(const Config& i_cfg, const std::string& i_key,
		            const std::string& i_section, const IOUtils::ThrowOptions& i_opt)
		            : proxycfg(i_cfg), key(i_key),section(i_section), opt(i_opt) { }

		template<typename T> operator T() {
			T tmp;
			proxycfg.getValue(key, section, tmp, opt);
			return tmp;
		}

		ConfigProxy& operator =(const ConfigProxy& /*i_cfg*/) {return *this;} //making VC++ happy...
};

} //end namespace mio

#endif
