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
#include <fstream>
#include <sstream>
#include <errno.h>
#include <cstring>

#include <meteoio/FileUtils.h>
#include <meteoio/meteoFilters/ProcessingBlock.h>
#include <meteoio/meteoFilters/FilterSuppr.h>
#include <meteoio/meteoFilters/FilterMin.h>
#include <meteoio/meteoFilters/FilterMax.h>
#include <meteoio/meteoFilters/FilterMinMax.h>
#include <meteoio/meteoFilters/FilterPotentialSW.h>
#include <meteoio/meteoFilters/FilterStdDev.h>
#include <meteoio/meteoFilters/FilterRate.h>
#include <meteoio/meteoFilters/FilterUnheatedPSUM.h>
#include <meteoio/meteoFilters/FilterTukey.h>
#include <meteoio/meteoFilters/FilterMAD.h>
#include <meteoio/meteoFilters/ProcAggregate.h>
#include <meteoio/meteoFilters/ProcIIR.h>
#include <meteoio/meteoFilters/ProcUndercatch_WMO.h>
#include <meteoio/meteoFilters/ProcUndercatch_Forland.h>
#include <meteoio/meteoFilters/ProcUndercatch_Hamon.h>
#include <meteoio/meteoFilters/ProcPSUMDistribute.h>
#include <meteoio/meteoFilters/ProcUnventilatedT.h>
#include <meteoio/meteoFilters/ProcShade.h>
#include <meteoio/meteoFilters/ProcAdd.h>
#include <meteoio/meteoFilters/ProcMult.h>
#include <meteoio/meteoFilters/ProcNoise.h>
#include <meteoio/meteoFilters/ProcExpSmoothing.h>
#include <meteoio/meteoFilters/ProcWMASmoothing.h>
#include <meteoio/meteoFilters/FilterNoChange.h>
#include <meteoio/meteoFilters/FilterTimeconsistency.h>
#include <meteoio/meteoFilters/FilterOffsetsnowdepth.h>
#include <meteoio/meteoFilters/FilterDeGrass.h>

namespace mio {
/**
 * @page processing Processing overview
 * The pre-processing infrastructure is described in ProcessingBlock (for its API). The goal of this page is to give an overview of the available filters and processing elements and their usage.
 *
 * @section processing_modes Modes of operation
 * It should be noted that filters often have two modes of operations: soft or hard. In soft mode, all value that is rejected is replaced by the filter parameter's value. This means that for a soft min filter set at 0.0, all values less than 0.0 will be replaced by 0.0. In hard mode, all rejected values are replaced by nodata.
 *
 * Several filter take arguments describing a processing window (for example, FilterStdDev). In such a case, two values are given:
 * - the minimum time span of the window
 * - the minimum number of points in the window
 *
 * The ProcessingBlock will walk through the data, starting at the current point and adding points to the processing window.  When both of these criterias are met,
 * the window is accepted. This means that a window defined as "6 21600" is a window that contains 6 points minimum AND spans at least 21600 seconds.
 *
 * @section processing_section Filtering section
 * The filters are specified for each parameter in the [Filters] section. This section contains
 * a list of the various meteo parameters (see MeteoData) with their associated choice of filtering algorithms and
 * optional parameters.The filters are applied serialy, in the order they are given in. An example of such section is given below:
 * @code
 * [Filters]
 * TA::filter1	= min_max
 * TA::arg1	= 230 330
 *
 * RH::filter1	= min_max
 * RH::arg1	= -0.2 1.2
 * RH::filter2	= min_max
 * RH::arg2	= soft 0.0 1.0
 *
 * PSUM::filter1	= min
 * PSUM::arg1	= -0.1
 * PSUM::filter2	= min
 * PSUM::arg2	= soft 0.
 * @endcode
 *
 * @section processing_available Available processing elements
 * New filters can easily be developed. The filters that are currently available are the following:
 * - NONE: this does nothing (this is useful in an \ref config_import "IMPORT" to overwritte previous filters);
 * - MIN: minimum check filter, see FilterMin
 * - MAX: maximum check filter, see FilterMax
 * - MIN_MAX: range check filter, see FilterMinMax
 * - RATE: rate of change filter, see FilterRate
 * - STD_DEV: reject data outside mean +/- k*stddev, see FilterStdDev
 * - MAD: median absolute deviation, see FilterMAD
 * - TUKEY: Tukey53H spike detection, based on median, see FilterTukey
 * - UNHEATED_RAINGAUGE: detection of snow melting in a rain gauge, see FilterUnheatedPSUM
 * - NO_CHANGE: reject data that changes too little (low variance), see FilterNoChange
 * - TIME_CONSISTENCY: reject data that changes too much , see FilterTimeconsistency
 * - DETECT_GRASS: detection of grass growing under the snow height sensor, see FilterDeGrass
 * - POTENTIALSW: ensuring physically realistic incoming short wave radiation, see FilterPotentialSW
 *
 * Some data transformations are also supported besides filtering, both very basic and generic data transformations:
 * - SUPPR: delete all or some data, see FilterSuppr
 * - ADD: adds a given offset to the data, see ProcAdd
 * - MULT: multiply the data by a given factor, see ProcMult
 * - NOISE: add to (or multiply by) randomly distributed noise, see ProcNoise
 *
 * As well as more specific data transformations:
 * - EXP_SMOOTHING: exponential smoothing of data, see ProcExpSmoothing
 * - WMA_SMOOTHING: weighted moving average smoothing of data, see ProcWMASmoothing
 * - IIR: Low Pass or High Pass critically damped filter, see ProcIIR
 * - AGGREGATE: various data aggregation algorithms, see ProcAggregate
 * - UNDERCATCH_WMO: WMO rain gauge correction for undercatch, using various correction models, see ProcUndercatch_WMO
 * - UNDERCATCH_FORLAND: Forland1996 rain gauge correction for solid and liquid undercatch, using various correction models, see ProcUndercatch_Forland
 * - UNDERCATCH_HAMON: Hamon1973 rain gauge correction for undercatch, see ProcUndercatch_Hamon
 * - UNVENTILATED_T: unventilated temperature sensor correction, see ProcUnventilatedT
 * - PSUM_DISTRIBUTE: distribute accumulated precipitation over preceeding timesteps, see ProcPSUMDistribute
 * - SHADE: apply a shading mask to the Incoming or Reflected Short Wave Radiation, see ProcShade
 *
 */

ProcessingBlock* BlockFactory::getBlock(const std::string& blockname, const std::vector<std::string>& vec_args, const Config& cfg)
{
	//the indenting is a little weird, this is in order to show the same groups as in the documentation above
	
	//normal filters
	if (blockname == "MIN"){
		return new FilterMin(vec_args, blockname);
	} else if (blockname == "MAX"){
		return new FilterMax(vec_args, blockname);
	} else if (blockname == "MIN_MAX"){
		return new FilterMinMax(vec_args, blockname);
	} else if (blockname == "RATE"){
		return new FilterRate(vec_args, blockname);
	} else if (blockname == "STD_DEV"){
		return new FilterStdDev(vec_args, blockname);
	} else if (blockname == "MAD"){
		return new FilterMAD(vec_args, blockname);
	} else if (blockname == "TUKEY"){
		return new FilterTukey(vec_args, blockname);
	} else if (blockname == "UNHEATED_RAINGAUGE"){
		return new FilterUnheatedPSUM(vec_args, blockname);
	} else if (blockname == "NO_CHANGE"){
		return new FilterNoChange(vec_args, blockname);
	} else if (blockname == "TIME_CONSISTENCY"){
		return new FilterTimeconsistency(vec_args, blockname);
	} else if (blockname == "OFFSNOW"){
		return new FilterOffsetsnowdepth(vec_args, blockname);
	} else if (blockname == "DETECT_GRASS"){
		return new FilterDeGrass(vec_args, blockname);
	} else if (blockname == "POTENTIALSW"){
		return new FilterPotentialSW(vec_args, blockname);
	}
	
	//general data transformations
	else if (blockname == "SUPPR"){
		return new FilterSuppr(vec_args, blockname, cfg.getConfigRootDir(), cfg.get("TIME_ZONE", "Input"));
	} else if (blockname == "ADD"){
		return new ProcAdd(vec_args, blockname, cfg.getConfigRootDir());
	} else if (blockname == "MULT"){
		return new ProcMult(vec_args, blockname, cfg.getConfigRootDir());
	} else if (blockname == "NOISE"){
		return new ProcNoise(vec_args, blockname);
	}
	
	//more specific data transformations
	else if (blockname == "EXP_SMOOTHING"){
		return new ProcExpSmoothing(vec_args, blockname);
	} else if (blockname == "WMA_SMOOTHING"){
		return new ProcWMASmoothing(vec_args, blockname);
	} else if (blockname == "IIR"){
		return new ProcIIR(vec_args, blockname);
	} else if (blockname == "AGGREGATE"){
		return new ProcAggregate(vec_args, blockname);
	} else if (blockname == "UNDERCATCH_WMO"){
		return new ProcUndercatch_WMO(vec_args, blockname);
	} else if (blockname == "UNDERCATCH_FORLAND"){
		return new ProcUndercatch_Forland(vec_args, blockname);
	} else if (blockname == "UNDERCATCH_HAMON"){
		return new ProcUndercatch_Hamon(vec_args, blockname);
	} else if (blockname == "UNVENTILATED_T"){
		return new ProcUnventilatedT(vec_args, blockname);
	} else if (blockname == "PSUM_DISTRIBUTE"){
		return new ProcPSUMDistribute(vec_args, blockname);
	} else if (blockname == "SHADE"){
		return new ProcShade(vec_args, blockname, cfg);
	} else {
		throw IOException("The processing block '"+blockname+"' does not exist! " , AT);
	}

}

const double ProcessingBlock::soil_albedo = .23; //grass
const double ProcessingBlock::snow_albedo = .85; //snow
const double ProcessingBlock::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

ProcessingBlock::ProcessingBlock(const std::string& name) : properties(), block_name(name)
{}

void ProcessingBlock::convert_args(const size_t& min_nargs, const size_t& max_nargs,
                               const std::vector<std::string>& vec_args, std::vector<double>& dbl_args) const
{
	const size_t nr_of_args = vec_args.size();
	if ((nr_of_args < min_nargs) || (nr_of_args > max_nargs))
		throw InvalidArgumentException("Wrong number of arguments for filter/processing element \"" + getName() + "\"", AT);

	for (size_t ii=0; ii<nr_of_args; ii++){
		double tmp = IOUtils::nodata;
		IOUtils::convertString(tmp, vec_args[ii]);
		dbl_args.push_back(tmp);
	}
}

bool ProcessingBlock::is_soft(std::vector<std::string>& vec_args) {
	for (size_t ii=0; ii<vec_args.size(); ++ii) {
		if (IOUtils::strToUpper( vec_args[ii] ) == "SOFT") {
			vec_args.erase( vec_args.begin() + ii );
			return true;
		}
	}

	return false;
}

std::string ProcessingBlock::getName() const {
	return block_name;
}

void ProcessingBlock::readCorrections(const std::string& filter, const std::string& filename, const char& c_type, const double& init, std::vector<double> &corrections)
{
	std::ifstream fin( filename.c_str() );
	if (fin.fail()) {
		std::ostringstream ss;
		ss << "Filter " << filter << ": ";
		ss << "error opening file \"" << filename << "\", possible reason: " << std::strerror(errno);
		throw AccessException(ss.str(), AT);
	}

	if (c_type=='m') corrections.resize(12, init);
	else if (c_type=='d') corrections.resize(366, init);
	else if (c_type=='h') corrections.resize(24, init);
	const size_t maxIndex = corrections.size();
	const size_t minIndex = (c_type=='h')? 0 : 1;

	const char eoln = FileUtils::getEoln(fin); //get the end of line character for the file

	try {
		size_t index, lcount=0;
		double value;
		do {
			lcount++;
			std::string line;
			getline(fin, line, eoln); //read complete line
			IOUtils::stripComments(line);
			IOUtils::trim(line);
			if (line.empty()) continue;

			std::istringstream iss( line );
			iss.setf(std::ios::fixed);
			iss.precision(std::numeric_limits<double>::digits10);
			iss >> std::skipws >> index;
			if ( !iss || index<minIndex || (index-minIndex)>=maxIndex) {
				std::ostringstream ss;
				ss << "Invalid index in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}
			iss >> std::skipws >> value;
			if ( !iss ){
				std::ostringstream ss;
				ss << "Invalid value in file " << filename << " at line " << lcount;
				throw InvalidArgumentException(ss.str(), AT);
			}

			corrections.at( index-minIndex ) = value;
		} while (!fin.eof());
		fin.close();
	} catch (const std::exception&){
		if (fin.is_open()) {//close fin if open
			fin.close();
		}
		throw;
	}
}

ProcessingBlock::~ProcessingBlock() {}

const std::string ProcessingBlock::toString() const {
	std::ostringstream os;
	os << "[" << block_name << " ";
	os << properties.toString();
	os << "]";
	return os.str();
}

const ProcessingProperties& ProcessingBlock::getProperties() const {
	return properties;
}

const std::string ProcessingProperties::toString() const
{
	std::ostringstream os;
	const double h_before = time_before.getJulian()*24.;
	const double h_after = time_after.getJulian()*24.;
	const size_t p_before = points_before;
	const size_t p_after = points_after;

	os << "{";
	if(h_before>0. || h_after>0.) os << "-" << h_before << " +" << h_after << " h; ";
	if(p_before>0 || p_after>0) os << "-" << p_before << " +" << p_after << " pts; ";
	if(stage==ProcessingProperties::first)
		os << "p¹";
	if(stage==ProcessingProperties::second)
		os << "p²";
	if(stage==ProcessingProperties::both)
		os << "p½";
	os << "}";
	return os.str();
}

} //end namespace
