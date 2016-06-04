/***********************************************************************************/
/*  Copyright 2012 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
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
#include <meteoio/meteoFilters/ProcMult.h>
#include <meteoio/FileUtils.h>

using namespace std;

namespace mio {

ProcMult::ProcMult(const std::vector<std::string>& vec_args, const std::string& name, const std::string& i_root_path)
         : ProcessingBlock(name), vecFactors(), root_path(i_root_path), factor(0.), type('c')
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcMult::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;

	if (type=='c') {
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			tmp *= factor;
		}
	} else if (type=='m') {
		int year, month, day;
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			ovec[ii].date.getDate(year, month, day);
			tmp *= vecFactors[ month-1 ]; //indices start at 0
		}
	} else if (type=='d') {
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			tmp *= vecFactors[ ovec[ii].date.getJulianDayNumber() ];
		}
	} else if (type=='h') {
		int year, month, day, hour;
		for (size_t ii=0; ii<ovec.size(); ii++){
			double& tmp = ovec[ii](param);
			if (tmp == IOUtils::nodata) continue; //preserve nodata values

			ovec[ii].date.getDate(year, month, day, hour);
			tmp *= vecFactors[ hour ];
		}
	}
}


void ProcMult::parse_args(const std::vector<std::string>& vec_args)
{
	const size_t nrArgs = vec_args.size();
	if (nrArgs==1) {
		type='c'; //constant
		if (!IOUtils::convertString(factor, vec_args[0]))
			throw InvalidArgumentException("Invalid factor \""+vec_args[0]+"\" specified for the "+getName()+" filter. If correcting for a period, please specify the period!", AT);
	} else if (nrArgs==2) {
		const string type_str=IOUtils::strToUpper( vec_args[0] );
		if (type_str=="MONTHLY") type='m';
		else if (type_str=="DAILY") type='d';
		else if (type_str=="HOURLY") type='h';
		else
			throw InvalidArgumentException("Invalid period \""+type_str+"\" specified for the "+getName()+" filter", AT);

		//if this is a relative path, prefix the path with the current path
		const std::string in_filename = vec_args[1];
		const std::string prefix = ( IOUtils::isAbsolutePath(in_filename) )? "" : root_path+"/";
		const std::string path = IOUtils::getPath(prefix+in_filename, true);  //clean & resolve path
		const std::string filename = path + "/" + IOUtils::getFilename(in_filename);
		ProcessingBlock::readCorrections(getName(), filename, type, 1., vecFactors);
	} else
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);
}

} //end namespace
