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
#include <ctime>
#include <cstdlib>

#include <meteoio/meteoFilters/ProcNoise.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteoLaws/Meteoconst.h>

using namespace std;

namespace mio {

ProcNoise::ProcNoise(const std::vector<std::string>& vec_args, const std::string& name)
          : ProcessingBlock(name), range(IOUtils::nodata), distribution(), type()
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first;
}

void ProcNoise::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	srand( static_cast<unsigned int>(time(NULL)) );
	ovec = ivec;
	
	if (type=='a' && distribution=='u') uniform_add(param, ovec);
	else if (type=='m' && distribution=='u') uniform_mult(param, ovec);
	else if (type=='a' && distribution=='n') normal_add(param, ovec);
	else if (type=='m' && distribution=='n') normal_mult(param, ovec);
}

 //add a number between -range/2 and +range/2
void ProcNoise::uniform_add(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	for (size_t ii=0; ii<ovec.size(); ii++){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		tmp += (2.*static_cast<double>(rand())/(RAND_MAX) - 1.) * range;
	}
}

 //add a number between -range*val/2 and +range*val/2
void ProcNoise::uniform_mult(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	for (size_t ii=0; ii<ovec.size(); ii++){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		tmp *= (1. + (2.*static_cast<double>(rand())/(RAND_MAX)-1.) * range);
	}
}

void ProcNoise::normal_add(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	for (size_t ii=0; ii<ovec.size(); ii++){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		tmp += getBoxMuller() * range;
	}
}

void ProcNoise::normal_mult(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	for (size_t ii=0; ii<ovec.size(); ii++){
		double& tmp = ovec[ii](param);
		if (tmp == IOUtils::nodata) continue; //preserve nodata values

		tmp += (1. + getBoxMuller() * range);
	}
}

// Box–Muller method for normally distributed random numbers
// see https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
// This generate a normally distributed signal of mean=0 and std_dev=1
// For numerical reasons, the extremes will always be less than 7 std_dev from the mean
double ProcNoise::getBoxMuller() const
{
	const double U = static_cast<double>(rand())/(RAND_MAX);
	const double V = static_cast<double>(rand())/(RAND_MAX);
	
	return  Optim::fastSqrt_Q3(-2.*log(U)) * cos(2.*Cst::PI*V);
}

void ProcNoise::parse_args(std::vector<std::string> vec_args)
{
	if (vec_args.size()!=3){
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName() + ": it requires a type, a distribution and a range", AT);
	}
	
	const string type_str=IOUtils::strToUpper( vec_args[0] );
	if (type_str=="ADD") 
		type='a';
	else if (type_str=="MULT") 
		type='m';
	else
		throw InvalidArgumentException("Invalid type \""+type_str+"\" specified for the "+getName()+" filter", AT);
		
	const string distribution_str=IOUtils::strToUpper( vec_args[1] );
	if (distribution_str=="UNIFORM")
		distribution='u';
	else if (distribution_str=="NORMAL")
		distribution='n';
	else
		throw InvalidArgumentException("Invalid distribution \""+distribution_str+"\" specified for the "+getName()+" filter", AT);
	
	if (!IOUtils::convertString(range, vec_args[2]))
		throw InvalidArgumentException("Invalid range specified for the "+getName()+" filter", AT);
}

} //end namespace
