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
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <errno.h>

#include <meteoio/meteoFilters/FilterSuppr.h>
#include <meteoio/FileUtils.h>

#include <fstream>

using namespace std;

namespace mio {

FilterSuppr::FilterSuppr(const std::vector<std::string>& vec_args, const std::string& name, const std::string& i_root_path, const double& i_TZ)
          : FilterBlock(name), suppr_dates(), root_path(i_root_path), TZ(i_TZ), range(IOUtils::nodata)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void FilterSuppr::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;
	if (ovec.empty()) return;
	
	if (!suppr_dates.empty()) {
		const std::string station_ID = ivec[0].meta.stationID;
		const std::map< std::string, std::set<Date> >::const_iterator station_it = suppr_dates.find( station_ID );
		if (station_it==suppr_dates.end()) return;
		
		for (size_t ii=0; ii<ovec.size(); ii++){
			const std::set<Date>::const_iterator it = station_it->second.find( ovec[ii].date );
			if (it!=station_it->second.end()) {
				ovec[ii](param) = IOUtils::nodata;
			}
		}
	} else if (range==IOUtils::nodata) { //remove all
		for (size_t ii=0; ii<ovec.size(); ii++){
			ovec[ii](param) = IOUtils::nodata;
		}
	} else { //only remove a given fraction
		const size_t set_size = ovec.size();
		const size_t nrRemove = static_cast<size_t>( round( (double)set_size*range ) );

		srand( static_cast<unsigned int>(time(NULL)) );
		size_t ii=1;
		while (ii<nrRemove) {
			const size_t idx = rand() % set_size;
			if (ivec[idx](param)!=IOUtils::nodata && ovec[idx](param)==IOUtils::nodata) continue; //the point was already removed
			
			ovec[idx](param)=IOUtils::nodata; //ie nodata points remain and are counted
			ii++;
		}
	}
}

void FilterSuppr::fillSuppr_dates(const std::string& filename)
{
	if (!FileUtils::validFileAndPath(filename)) throw InvalidNameException(filename, AT);
	if (!FileUtils::fileExists(filename)) throw NotFoundException(filename, AT);
	
	std::ifstream fin(filename.c_str());
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << block_name << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}
	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file
	
	try {
		size_t lcount=0;
		do {
			lcount++;
			std::string line;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;
			
			std::vector<std::string> vecString;
			const size_t nrElems = IOUtils::readLineToVec(line, vecString);
			if (nrElems!=2) {
				std::ostringstream os;
				os << "Invalid syntax for filter " << block_name << " in file \"" << filename << "\": expecting 2 arguments, got " << nrElems;
				throw InvalidFormatException(os.str(), AT);
			}
			
			Date d1;
			if (!IOUtils::convertString(d1, vecString[1], TZ))
				throw InvalidFormatException("Could not process date "+vecString[1]+" in file \""+filename+"\"", AT);
			
			const std::string station_ID = vecString[0];
			suppr_dates[ station_ID ].insert( d1 );
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
}

void FilterSuppr::parse_args(std::vector<std::string> vec_args) 
{
	const size_t nrArgs = vec_args.size();

	if (nrArgs>1)
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);
	
	if (nrArgs==1) {
		if (IOUtils::isNumeric(vec_args[0])) {
			if (!IOUtils::convertString(range, vec_args[0]))
				throw InvalidArgumentException("Invalid range \""+vec_args[0]+"\" specified for the "+getName()+" filter.", AT);
			if (range<0. || range>1.)
				throw InvalidArgumentException("Wrong range for filter " + getName() + ", it should be between 0 and 1", AT);
		} else {
			const std::string in_filename( vec_args[0] );
			const std::string prefix = ( FileUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
			const std::string path = FileUtils::getPath(prefix+in_filename, true);  //clean & resolve path
			const std::string filename = path + "/" + FileUtils::getFilename(in_filename);
		
			fillSuppr_dates(filename);
		}
	}
}

} //end namespace
