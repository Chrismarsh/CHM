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
#include <meteoio/dataClasses/MeteoData.h>
#include <meteoio/dataClasses/StationData.h>

#include <cmath>
#include <limits>

using namespace std;
namespace mio {

/************************************************************
 * static section                                           *
 ************************************************************/
const size_t MeteoGrids::nrOfParameters =  MeteoGrids::lastparam - MeteoGrids::firstparam + 1;
std::vector<std::string> MeteoGrids::paramname;
const bool MeteoGrids::__init = MeteoGrids::initStaticData();

bool MeteoGrids::initStaticData()
{
	//the order must be the same as in the enum
	paramname.push_back("TA");
	paramname.push_back("RH");
	paramname.push_back("QI");
	paramname.push_back("TD");
	paramname.push_back("VW");
	paramname.push_back("DW");
	paramname.push_back("VW_MAX");
	paramname.push_back("ISWR");
	paramname.push_back("RSWR");
	paramname.push_back("ISWR_DIFF");
	paramname.push_back("ISWR_DIR");
	paramname.push_back("ILWR");
	paramname.push_back("TAU_CLD");
	paramname.push_back("HS");
	paramname.push_back("PSUM");
	paramname.push_back("PSUM_PH");
	paramname.push_back("PSUM_L");
	paramname.push_back("PSUM_S");
	paramname.push_back("TSG");
	paramname.push_back("TSS");
	paramname.push_back("P");
	paramname.push_back("P_SEA");
	paramname.push_back("U");
	paramname.push_back("V");
	paramname.push_back("W");
	paramname.push_back("SWE");
	paramname.push_back("ROT");
	paramname.push_back("ALB");
	paramname.push_back("DEM");
	paramname.push_back("SHADE");
	paramname.push_back("SLOPE");
	paramname.push_back("AZI");

	return true;
}

const std::string& MeteoGrids::getParameterName(const size_t& parindex)
{
	if (parindex >= MeteoGrids::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to get name for parameter that does not exist", AT);

	return paramname[parindex];
}

/************************************************************
 * static section                                           *
 ************************************************************/
const double MeteoData::epsilon = 1e-5;
const size_t MeteoData::nrOfParameters =  MeteoData::lastparam - MeteoData::firstparam + 1;
map<size_t, string> MeteoData::static_meteoparamname;
std::vector<std::string> MeteoData::s_default_paramname;
const bool MeteoData::__init = MeteoData::initStaticData();

bool MeteoData::initStaticData()
{
	//Associate unsigned int value and a string representation of a meteo parameter
	static_meteoparamname[P]      = "P";
	static_meteoparamname[TA]     = "TA";
	static_meteoparamname[RH]     = "RH";
	static_meteoparamname[TSG]    = "TSG";
	static_meteoparamname[TSS]    = "TSS";
	static_meteoparamname[HS]     = "HS";
	static_meteoparamname[VW]     = "VW";
	static_meteoparamname[DW]     = "DW";
	static_meteoparamname[VW_MAX] = "VW_MAX";
	static_meteoparamname[RSWR]   = "RSWR";
	static_meteoparamname[ISWR]   = "ISWR";
	static_meteoparamname[ILWR]   = "ILWR";
	static_meteoparamname[TAU_CLD]= "TAU_CLD";
	static_meteoparamname[PSUM]    = "PSUM";
	static_meteoparamname[PSUM_PH]    = "PSUM_PH";

	s_default_paramname.push_back("P");
	s_default_paramname.push_back("TA");
	s_default_paramname.push_back("RH");
	s_default_paramname.push_back("TSG");
	s_default_paramname.push_back("TSS");
	s_default_paramname.push_back("HS");
	s_default_paramname.push_back("VW");
	s_default_paramname.push_back("DW");
	s_default_paramname.push_back("VW_MAX");
	s_default_paramname.push_back("RSWR");
	s_default_paramname.push_back("ISWR");
	s_default_paramname.push_back("ILWR");
	s_default_paramname.push_back("TAU_CLD");
	s_default_paramname.push_back("PSUM");
	s_default_paramname.push_back("PSUM_PH");

	return true;
}

const std::string& MeteoData::getParameterName(const size_t& parindex)
{
	if (parindex >= MeteoData::nrOfParameters)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);

	return MeteoData::static_meteoparamname[parindex];
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

const std::string& MeteoData::getNameForParameter(const size_t& parindex) const
{
	if (parindex >= nrOfAllParameters)
		throw IndexOutOfBoundsException("Trying to get name for parameter that does not exist", AT);

	return param_name[parindex];
}

bool MeteoData::param_exists(const std::string& i_paramname) const
{
	const size_t current_size = param_name.size();
	for (size_t ii = 0; ii<current_size; ii++) {
		if (param_name[ii] == i_paramname)
			return true;
	}

	return false;
}

size_t MeteoData::addParameter(const std::string& i_paramname)
{
	//check if name is already taken
	const size_t current_index = getParameterIndex(i_paramname);
	if (current_index != IOUtils::npos)
		return current_index; //do nothing, because parameter is already present

	//add parameter
	param_name.push_back(i_paramname);
	data.push_back(IOUtils::nodata);

	//Increase nrOfAllParameters
	nrOfAllParameters++;

	return (nrOfAllParameters - 1);
}

size_t MeteoData::getNrOfParameters() const
{
	return nrOfAllParameters;
}

MeteoData::MeteoData()
         : date(0.0, 0.), meta(), param_name(s_default_paramname), data(MeteoData::nrOfParameters, IOUtils::nodata), nrOfAllParameters(MeteoData::nrOfParameters), resampled(false)
{ }

MeteoData::MeteoData(const Date& date_in)
         : date(date_in), meta(), param_name(s_default_paramname), data(MeteoData::nrOfParameters, IOUtils::nodata), nrOfAllParameters(MeteoData::nrOfParameters), resampled(false)
{ }

MeteoData::MeteoData(const Date& date_in, const StationData& meta_in)
         : date(date_in), meta(meta_in), param_name(s_default_paramname), data(MeteoData::nrOfParameters, IOUtils::nodata), nrOfAllParameters(MeteoData::nrOfParameters), resampled(false)
{ }

void MeteoData::setDate(const Date& in_date)
{
	date = in_date;
}

void MeteoData::reset()
{
	std::fill(data.begin(), data.end(), IOUtils::nodata);
}

/**
* @brief Standardize nodata values
* The plugin-specific nodata values are replaced by MeteoIO's internal nodata value
*/
void MeteoData::standardizeNodata(const double& plugin_nodata) {
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		//loop through all meteo params and check whether they're nodata values
		if (data[ii] <= plugin_nodata)
			data[ii] = IOUtils::nodata;
	}
}

bool MeteoData::isResampled() const
{
	return resampled;
}

void MeteoData::setResampled(const bool& in_resampled)
{
	resampled = in_resampled;
}

bool MeteoData::operator==(const MeteoData& in) const
{
	//An object is equal if the date is equal and all meteo parameters are equal
	if (date != in.date) {
		return false;
	}

	if (nrOfAllParameters != in.nrOfAllParameters) { //the number of meteo parameters has to be consistent
		return false;
	}

	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		//const double epsilon_rel = (fabs(data[ii]) < fabs(in.data[ii]) ? fabs(in.data[ii]) : fabs(data[ii])) * std::numeric_limits<double>::epsilon(); // Hack not working...
		//const double epsilon_rel = (fabs(data[ii]) < fabs(in.data[ii]) ? fabs(in.data[ii]) : fabs(data[ii])) * 0.0000001; // Hack not working with 0 == 0 ....
		if( !IOUtils::checkEpsilonEquality(data[ii], in.data[ii], epsilon) ){
			return false;
		}
	}

	return true;
}

bool MeteoData::operator!=(const MeteoData& in) const
{
	return !(*this==in);
}

double& MeteoData::operator()(const size_t& parindex)
{
#ifndef NOSAFECHECKS
	if (parindex >= nrOfAllParameters)//getNrOfParameters())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return data[parindex];
}

const double& MeteoData::operator()(const size_t& parindex) const
{
#ifndef NOSAFECHECKS
	if (parindex >= nrOfAllParameters)//getNrOfParameters())
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist", AT);
#endif
	return data[parindex];
}

double& MeteoData::operator()(const std::string& parname)
{
	const size_t index = getParameterIndex(parname);
	if (index == IOUtils::npos)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist: " + parname, AT);

	return operator()(index);
}

const double& MeteoData::operator()(const std::string& parname) const
{
	const size_t index = getParameterIndex(parname);
	if (index == IOUtils::npos)
		throw IndexOutOfBoundsException("Trying to access meteo parameter that does not exist: " + parname, AT);

	return operator()(index);
}

size_t MeteoData::getParameterIndex(const std::string& parname) const
{
	for (size_t ii=0; ii<nrOfAllParameters; ii++) {
		if (param_name[ii] == parname)
			return ii;
	}

	return IOUtils::npos; //parameter not a part of MeteoData
}

const std::string MeteoData::toString() const {
	std::ostringstream os;
	os << "<meteo>\n";
	os << meta.toString();
	os << date.toString(Date::FULL) << "\n";

	for (size_t ii=0; ii<getNrOfParameters(); ii++) {
		const double& value = operator()(ii);
		if (value != IOUtils::nodata)
			os << setw(8) << getNameForParameter(ii) << ":" << setw(15) << value << endl;
	}

	os << "</meteo>\n";
	return os.str();
}

std::iostream& operator<<(std::iostream& os, const MeteoData& data) {
	os << data.date;
	os << data.meta;
	const size_t s_vector = data.param_name.size();
	os.write(reinterpret_cast<const char*>(&s_vector), sizeof(size_t));
	for (size_t ii=0; ii<s_vector; ii++) {
		const size_t s_string = data.param_name[ii].size();
		os.write(reinterpret_cast<const char*>(&s_string), sizeof(size_t));
		os.write(reinterpret_cast<const char*>(&data.param_name[ii][0]), s_string*sizeof(data.param_name[ii][0]));
	}

	const size_t s_data = data.data.size();
	os.write(reinterpret_cast<const char*>(&s_data), sizeof(s_data));
	os.write(reinterpret_cast<const char*>(&data.data[0]), s_data*sizeof(data.data[0]));

	os.write(reinterpret_cast<const char*>(&data.nrOfAllParameters), sizeof(data.nrOfAllParameters));
	os.write(reinterpret_cast<const char*>(&data.resampled), sizeof(data.resampled));
	return os;
}

std::iostream& operator>>(std::iostream& is, MeteoData& data) {
	is >> data.date;
	is >> data.meta;
	size_t s_vector;
	is.read(reinterpret_cast<char*>(&s_vector), sizeof(size_t));
	data.param_name.resize(s_vector);
	for (size_t ii=0; ii<s_vector; ii++) {
		size_t s_string;
		is.read(reinterpret_cast<char*>(&s_string), sizeof(size_t));
		data.param_name[ii].resize(s_string);
		is.read(reinterpret_cast<char*>(&data.param_name[ii][0]), s_string*sizeof(data.param_name[ii][0]));
	}

	size_t s_data;
	is.read(reinterpret_cast<char*>(&s_data), sizeof(size_t));
	data.data.resize(s_data);
	is.read(reinterpret_cast<char*>(&data.data[0]), s_data*sizeof(data.data[0]));

	is.read(reinterpret_cast<char*>(&data.nrOfAllParameters), sizeof(data.nrOfAllParameters));
	is.read(reinterpret_cast<char*>(&data.resampled), sizeof(data.resampled));
	return is;
}

void MeteoData::merge(std::vector<MeteoData>& vec1, const std::vector<MeteoData>& vec2, const bool& simple_merge)
{
	if (vec2.empty()) return;

	if (simple_merge || vec1.empty()) {
		vec1.reserve( vec1.size()+vec2.size() );
		for(size_t ii=0; ii<vec2.size(); ii++) vec1.push_back( vec2[ii] );
	} else {
		for (size_t ii=0; ii<vec2.size(); ii++) merge(vec1, vec2[ii]);
	}
}

void MeteoData::merge(std::vector<MeteoData>& vec, const MeteoData& meteo2, const bool& simple_merge)
{
	if (!simple_merge) {
		for (size_t ii=0; ii<vec.size(); ii++) {
			//two stations are considered the same if they point to the same 3D position
			if (vec[ii].meta.position==meteo2.meta.position) {
				vec[ii].merge(meteo2);
				return;
			}
		}
	}

	//the station was not found in the vector -> adding it
	vec.push_back( meteo2 );
}

void MeteoData::merge(std::vector<MeteoData>& vec)
{
	const size_t nElems = vec.size();
	if (nElems<2) return;
	
	std::vector<MeteoData> vecResult;
	std::vector<size_t> mergeIdx(nElems);
	for(size_t ii=0; ii<nElems; ii++) mergeIdx[ii] = ii;
	
	for (size_t ii=0; ii<nElems; ii++) {
		if (mergeIdx[ii]==IOUtils::npos) continue; //this element has already been merged, skip
		for (size_t jj=ii+1; jj<nElems; jj++) {
			if (vec[ii].meta.position==vec[jj].meta.position) {
				vec[ii].merge( vec[jj] );
				mergeIdx[jj]=IOUtils::npos; //this element will be skipped in the next loops
			}
		}
		vecResult.push_back( vec[ii] );
	}
	
	vec.swap( vecResult );
}

MeteoData MeteoData::merge(const MeteoData& meteo1, const MeteoData& meteo2)
{
	MeteoData tmp(meteo1);
	tmp.merge(meteo2);
	return tmp;
}

void MeteoData::merge(const MeteoData& meteo2)
{
	if (!date.isUndef() && !meteo2.date.isUndef() && date!=meteo2.date) {
		//the data must be time synchronized!
		std::ostringstream ss;
		ss << "Trying to merge MeteoData at " << date.toString(Date::ISO);
		ss << " with MeteoData at " << meteo2.date.toString(Date::ISO);
		throw InvalidArgumentException(ss.str(), AT);
	}

	if (date.isUndef()) date=meteo2.date;
	meta.merge(meteo2.meta);

	if (meteo2.resampled==true ) resampled=true;

	//merge standard parameters
	for (size_t ii=0; ii<nrOfParameters; ii++) {
		if (data[ii]==IOUtils::nodata) data[ii]=meteo2.data[ii];
	}

	//merge extra parameters
	const size_t nrExtra1 = nrOfAllParameters - nrOfParameters;
	const size_t nrExtra2 = meteo2.nrOfAllParameters - nrOfParameters;

	//no extra parameters to add -> return
	if (nrExtra2==0) return;

	//extra parameters only in meteo2 -> add them
	if (nrExtra1==0) {
		for (size_t ii=0; ii<nrExtra2; ii++) {
			const size_t new_idx = addParameter( meteo2.param_name[nrOfParameters+ii] );
			data[new_idx] = meteo2.data[nrOfParameters+ii];
		}
		return;
	}

	//extra parameters in both the current meteodata and meteo2 -> tedious merge...
	std::vector<bool> meteo2_flags(nrExtra2, false); //to keep track of which elements have been copied
	//merge the extra field of the current meteodata
	for (size_t ii=0; ii<nrExtra1; ii++) {
		const std::string curr_name = param_name[nrOfParameters+ii];
		//look for this parameter in meteo2
		for (size_t jj=0; jj<nrExtra2; jj++) {
			if (meteo2.param_name[nrOfParameters+jj]==curr_name ) {
				meteo2_flags[jj] = true;
				if (data[nrOfParameters+ii]==IOUtils::nodata) data[nrOfParameters+ii]=meteo2.data[nrOfParameters+jj];
				break;
			}
		}
	}
	//merge the extra fields of meteo2 that were NOT in the current meteodata
	for (size_t ii=0; ii<nrExtra2; ii++) {
		if (meteo2_flags[ii]==false) {
			const size_t new_idx = addParameter( meteo2.param_name[nrOfParameters+ii] );
			data[new_idx] = meteo2.data[nrOfParameters+ii];
		}
	}
}

} //namespace
