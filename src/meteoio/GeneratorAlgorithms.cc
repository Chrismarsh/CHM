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
#include <meteoio/GeneratorAlgorithms.h>
#include <meteoio/MathOptim.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <meteoio/meteoLaws/Meteoconst.h>
#include <meteoio/meteoFilters/ProcPSUMDistribute.h> //for the precipitation distribution

#include <algorithm>

using namespace std;

namespace mio {

GeneratorAlgorithm* GeneratorAlgorithmFactory::getAlgorithm(const std::string& i_algoname, const std::vector<std::string>& vecArgs)
{
	std::string algoname(i_algoname);
	IOUtils::toUpper(algoname);

	if (algoname == "CST"){
		return new ConstGenerator(vecArgs, i_algoname);
	} else if (algoname == "SIN"){
		return new SinGenerator(vecArgs, i_algoname);
	} else if (algoname == "STD_PRESS"){
		return new StandardPressureGenerator(vecArgs, i_algoname);
	} else if (algoname == "RELHUM"){
		return new RhGenerator(vecArgs, i_algoname);
	} else if (algoname == "TS_OLWR"){
		return new TsGenerator(vecArgs, i_algoname);
	} else if (algoname == "ISWR_ALBEDO"){
		return new IswrAlbedoGenerator(vecArgs, i_algoname);
	} else if (algoname == "CLEARSKY_LW"){
		return new ClearSkyLWGenerator(vecArgs, i_algoname);
	} else if (algoname == "ALLSKY_LW"){
		return new AllSkyLWGenerator(vecArgs, i_algoname);
	} else if (algoname == "ALLSKY_SW"){
		return new AllSkySWGenerator(vecArgs, i_algoname);
	} else if (algoname == "ESOLIP"){
		return new ESOLIPGenerator(vecArgs, i_algoname);
	} else if (algoname == "PPHASE"){
		return new PPhaseGenerator(vecArgs, i_algoname);
	} else {
		throw IOException("The generator algorithm '"+algoname+"' is not implemented" , AT);
	}
}

std::string GeneratorAlgorithm::getAlgo() const {
	return algo;
}

void GeneratorAlgorithm::parse_args(const std::vector<std::string>& vecArgs)
{
	if(!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

////////////////////////////////////////////////////////////////////////

void ConstGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(vecArgs.size()==1) {
		IOUtils::convertString(constant, vecArgs[0]);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool ConstGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata)
		value = constant;

	return true; //all missing values could be filled
}

bool ConstGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		generate(param, vecMeteo[ii]);
	}

	return true; //all missing values could be filled
}

void SinGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(vecArgs.size()==4) {
		const string type_str=IOUtils::strToUpper(vecArgs[0]);
		if( type_str=="YEARLY" ) type='y';
		else if( type_str=="DAILY" ) type='d';
		else
			throw InvalidArgumentException("Invalid period \""+type_str+"\" specified for the "+algo+" generator", AT);

		double min, max;
		IOUtils::convertString(min, vecArgs[1]);
		IOUtils::convertString(max, vecArgs[2]);
		amplitude = 0.5*(max-min); //the user provides min, max
		offset = min+amplitude;
		IOUtils::convertString(phase, vecArgs[3]);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool SinGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		double t; //also, the minimum must occur at 0 if phase=0
		if(type=='y') {
			t = (static_cast<double>(md.date.getJulianDayNumber()) - phase*365.25) / 366.25 - .25;
		} else if(type=='d') {
			const double julian = md.date.getJulian();
			t = (julian - Optim::intPart(julian) - phase) + .25; //watch out: julian day starts at noon!
		} else {
			std::ostringstream ss;
			ss << "Invalid period \"" << type << "\" specified for the " << algo << " generator";
			throw InvalidArgumentException(ss.str(), AT);
		}

		const double w = 2.*Cst::PI;
		value = amplitude * sin(w*t) + offset;
	}

	return true; //all missing values could be filled
}

bool SinGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		generate(param, vecMeteo[ii]);
	}

	return true; //all missing values could be filled
}


bool StandardPressureGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		const double altitude = md.meta.position.getAltitude();
		if(altitude==IOUtils::nodata) return false;
		value = Atmosphere::stdAirPressure(altitude);
	}

	return true; //all missing values could be filled
}

bool StandardPressureGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	const double altitude = vecMeteo.front().meta.position.getAltitude(); //if the stations move, this has to be in the loop
	if(altitude==IOUtils::nodata) return false;

	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if(value == IOUtils::nodata)
			value = Atmosphere::stdAirPressure(altitude);
	}

	return true; //all missing values could be filled
}


bool RhGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		const double TA = md(MeteoData::TA);
		if (TA==IOUtils::nodata) //nothing else we can do here
			return false;

		//first chance to compute RH
		if (md.param_exists("TD")) {
			const double TD = md("TD");
			if (TD!=IOUtils::nodata)
				value = Atmosphere::DewPointtoRh(TD, TA, false);
		}

		//second chance to try to compute RH
		if (value==IOUtils::nodata && md.param_exists("SH")) {
			const double SH = md("SH");
			const double altitude = md.meta.position.getAltitude();
			if (SH!=IOUtils::nodata && altitude!=IOUtils::nodata)
				value = Atmosphere::specToRelHumidity(altitude, TA, SH);
		}

		if (value==IOUtils::nodata) return false;
	}

	return true; //all missing values could be filled
}

bool RhGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	const double altitude = vecMeteo.front().meta.position.getAltitude(); //if the stations move, this has to be in the loop

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if(value == IOUtils::nodata) {
			const double TA = vecMeteo[ii](MeteoData::TA);
			if (TA==IOUtils::nodata) { //nothing else we can do here
				all_filled=false;
				continue;
			}

			//first chance to compute RH
			if (vecMeteo[ii].param_exists("TD")) {
				const double TD = vecMeteo[ii]("TD");
				if (TD!=IOUtils::nodata)
					value = Atmosphere::DewPointtoRh(TD, TA, false);
			}

			//second chance to try to compute RH
			if (value==IOUtils::nodata && vecMeteo[ii].param_exists("SH")) {
				const double SH = vecMeteo[ii]("SH");
				if (SH!=IOUtils::nodata && altitude!=IOUtils::nodata)
					value = Atmosphere::specToRelHumidity(altitude, TA, SH);
			}

			if (value==IOUtils::nodata) all_filled=false;
		}
	}

	return all_filled;
}


const double TsGenerator::e_snow = .983; //snow emissivity (0.969 - 0.997)
const double TsGenerator::e_soil = .9805; //grass emissivity (0.975 - 0.986)
const double TsGenerator::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

bool TsGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		const double olwr = md("OLWR");
		if (olwr==IOUtils::nodata) //nothing else we can do here
			return false;

		const double hs = md(MeteoData::HS);
		const double ea = (hs==IOUtils::nodata)? .5*(e_snow+e_soil) : (hs>snow_thresh)? e_snow : e_soil;

		//value = pow( olwr / ( ea * Cst::stefan_boltzmann ), 0.25);
		value = Optim::invSqrt( Optim::invSqrt(olwr / ( ea * Cst::stefan_boltzmann )) );

		if (value==IOUtils::nodata) return false;
	}

	return true; //all missing values could be filled
}

bool TsGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	bool status = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (!generate(param, vecMeteo[ii]))
			status = false;
	}

	return status;
}


const double IswrAlbedoGenerator::soil_albedo = .23; //grass
const double IswrAlbedoGenerator::snow_albedo = .85; //snow
const double IswrAlbedoGenerator::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

bool IswrAlbedoGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		const double HS=md(MeteoData::HS), RSWR=md(MeteoData::RSWR), ISWR=md(MeteoData::ISWR);

		double albedo = .5;
		if (HS!=IOUtils::nodata)
			albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;

		if (param==MeteoData::ISWR && (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)) {
			value = RSWR / albedo;
			return true;
		}

		if (param==MeteoData::RSWR && (ISWR!=IOUtils::nodata && HS!=IOUtils::nodata)) {
			value = ISWR * albedo;
			return true;
		}

		return false;
	}

	return true; //all missing values could be filled
}

bool IswrAlbedoGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	bool status = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (!generate(param, vecMeteo[ii]))
			status = false;
	}

	return status;
}


void ClearSkyLWGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(vecArgs.size()==1) {
		const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);

		if (user_algo=="BRUTSAERT") model = BRUTSAERT;
		else if (user_algo=="DILLEY") model = DILLEY;
		else if (user_algo=="PRATA") model = PRATA;
		else if (user_algo=="CLARK") model = CLARK;
		else if (user_algo=="TANG") model = TANG;
		else if (user_algo=="IDSO") model = IDSO;
		else
			throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for the "+algo+" generator", AT);
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool ClearSkyLWGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH);
		if (TA==IOUtils::nodata || RH==IOUtils::nodata) return false;

		if (model==BRUTSAERT)
			value = Atmosphere::Brutsaert_ilwr(RH, TA);
		else if (model==DILLEY)
			value = Atmosphere::Dilley_ilwr(RH, TA);
		else if (model==PRATA)
			value = Atmosphere::Prata_ilwr(RH, TA);
		else if (model==CLARK)
			value = Atmosphere::Clark_ilwr(RH, TA);
		else if (model==TANG)
			value = Atmosphere::Tang_ilwr(RH, TA);
		else if (model==IDSO)
			value = Atmosphere::Idso_ilwr(RH, TA);
	}

	return true; //all missing values could be filled
}

bool ClearSkyLWGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		const bool status = generate(param, vecMeteo[ii]);
		if(status==false) all_filled=false;
	}

	return all_filled;
}


const double AllSkyLWGenerator::soil_albedo = .23; //grass
const double AllSkyLWGenerator::snow_albedo = .85; //snow
const double AllSkyLWGenerator::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo

void AllSkyLWGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if (vecArgs.size()==1) {
		const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);

		if (user_algo=="OMSTEDT") model = OMSTEDT;
		else if (user_algo=="KONZELMANN") model = KONZELMANN;
		else if (user_algo=="UNSWORTH") model = UNSWORTH;
		else if (user_algo=="CRAWFORD") model = CRAWFORD;
		else
			throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied for the "+algo+" generator", AT);

		if (model==CRAWFORD) clf_model = CLF_CRAWFORD;
	} else { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

double AllSkyLWGenerator::getCloudiness(const MeteoData& md, SunObject& sun, bool &is_night)
{
	//we know that TA and RH are available, otherwise we would not get called
	const double TA=md(MeteoData::TA), RH=md(MeteoData::RH), HS=md(MeteoData::HS), RSWR=md(MeteoData::RSWR);
	double ISWR=md(MeteoData::ISWR);

	is_night = false;

	double albedo = .5;
	if (RSWR==IOUtils::nodata || ISWR==IOUtils::nodata || RSWR<=0 || ISWR<=0) {
		if (HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
			albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;

		if (ISWR==IOUtils::nodata && (RSWR!=IOUtils::nodata && HS!=IOUtils::nodata)) {
			ISWR = RSWR / albedo;
		}
	} else {
		albedo = RSWR / ISWR;
		if (albedo>=1.) albedo=0.99;
		if (albedo<=0.) albedo=0.01;
	}

	if (ISWR<5.) {
		is_night = true;
		return IOUtils::nodata;
	}

	if (ISWR==IOUtils::nodata) return IOUtils::nodata; //no way to get ISWR

	sun.calculateRadiation(TA, RH, albedo);
	double toa, direct, diffuse;
	sun.getHorizontalRadiation(toa, direct, diffuse);
	const double iswr_clear_sky = direct+diffuse;

	//at sunrise or sunset, we might get clf<0 or clf>1 -> return nodata in order to use interpolation instead
	if (iswr_clear_sky<5. || iswr_clear_sky<ISWR) {
		is_night = true;
		return IOUtils::nodata;
	}

	if (clf_model==KASTEN) {
		const double clf = Atmosphere::Kasten_cloudiness(ISWR/iswr_clear_sky);
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else if (clf_model==CLF_CRAWFORD) {
		const double clf = 1. - ISWR/iswr_clear_sky;
		if (clf<0. || clf>1.) return IOUtils::nodata;
		return clf;
	} else
		return IOUtils::nodata; //this should never happen
}

bool AllSkyLWGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA), RH=md(MeteoData::RH), TAU_CLD=md(MeteoData::TAU_CLD);
		if (TA==IOUtils::nodata || RH==IOUtils::nodata) return false;
		double cloudiness = (TAU_CLD!=IOUtils::nodata)? Atmosphere::Kasten_cloudiness( 1.-TAU_CLD ) : IOUtils::nodata;

		const string station_hash = md.meta.stationID + ":" + md.meta.stationName;
		const double julian_gmt = md.date.getJulian(true);
		bool cloudiness_from_cache = false;

		//try to get a cloudiness value
		if (cloudiness==IOUtils::nodata) {
			const double lat = md.meta.position.getLat();
			const double lon = md.meta.position.getLon();
			const double alt = md.meta.position.getAltitude();
			SunObject sun;
			sun.setLatLon(lat, lon, alt);
			sun.setDate(julian_gmt, 0.);

			bool is_night;
			cloudiness = getCloudiness(md, sun, is_night);
			if (cloudiness==IOUtils::nodata && !is_night) return false;

			if (is_night) { //interpolate the cloudiness over the night
				const map< string, pair<double, double> >::const_iterator it = last_cloudiness.find(station_hash);
				if (it==last_cloudiness.end()) return false;

				cloudiness_from_cache = true;
				const double last_cloudiness_julian = it->second.first;
				const double last_cloudiness_value = it->second.second;
				if ((julian_gmt - last_cloudiness_julian) < 1.) cloudiness = last_cloudiness_value;
				else return false;
			}
		}

		//run the ILWR parametrization
		if (model==OMSTEDT)
			value = Atmosphere::Omstedt_ilwr(RH, TA, cloudiness);
		else if (model==KONZELMANN)
			value = Atmosphere::Konzelmann_ilwr(RH, TA, cloudiness);
		else if (model==UNSWORTH)
			value = Atmosphere::Unsworth_ilwr(RH, TA, IOUtils::nodata, IOUtils::nodata, cloudiness);
		else if (model==CRAWFORD) {
			int year, month, day;
			md.date.getDate(year, month, day);
			value = Atmosphere::Crawford_ilwr(RH, TA, IOUtils::nodata, IOUtils::nodata, static_cast<unsigned char>(month), cloudiness);
		}

		//save the last valid cloudiness
		if (!cloudiness_from_cache)
			last_cloudiness[station_hash] = pair<double,double>( julian_gmt, cloudiness );
	}

	return true; //all missing values could be filled
}

bool AllSkyLWGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		const bool status = generate(param, vecMeteo[ii]);
		if(status==false) all_filled=false;
	}

	return all_filled;
}


const double AllSkySWGenerator::soil_albedo = .23; //grass
const double AllSkySWGenerator::snow_albedo = .85; //snow
const double AllSkySWGenerator::snow_thresh = .1; //if snow height greater than this threshold -> snow albedo
void AllSkySWGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	//Get the optional arguments for the algorithm: constant value to use
	if(!vecArgs.empty()) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool AllSkySWGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if(value == IOUtils::nodata) {
		const double ISWR=md(MeteoData::ISWR), RSWR=md(MeteoData::RSWR), HS=md(MeteoData::HS), TAU_CLD=md(MeteoData::TAU_CLD);
		double TA=md(MeteoData::TA), RH=md(MeteoData::RH), ILWR=md(MeteoData::ILWR);

		const double lat = md.meta.position.getLat();
		const double lon = md.meta.position.getLon();
		const double alt = md.meta.position.getAltitude();
		if(lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;

		double albedo = .5;
		if(RSWR==IOUtils::nodata || ISWR==IOUtils::nodata) {
			if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
				albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
		} else if(ISWR>0. && RSWR>0.) { //this could happen if the user calls this generator for a copy parameter, etc
			albedo = RSWR / ISWR;
			if(albedo>=1.) albedo=0.99;
			if(albedo<=0.) albedo=0.01;
		}

		if(TA==IOUtils::nodata || RH==IOUtils::nodata) {
			//set TA & RH so the reduced precipitable water will get an average value
			TA=274.98;
			RH=0.666;
			ILWR=IOUtils::nodata; //skip solarIndex correction
		}

		sun.setLatLon(lat, lon, alt);
		sun.setDate(md.date.getJulian(true), 0.);
		const double solarIndex = (TAU_CLD!=IOUtils::nodata)? TAU_CLD : (ILWR!=IOUtils::nodata)? getSolarIndex(TA, RH, ILWR) : 1.;

		const double P=md(MeteoData::P);
		if(P==IOUtils::nodata)
			sun.calculateRadiation(TA, RH, albedo);
		else
			sun.calculateRadiation(TA, RH, P, albedo);

		double toa, direct, diffuse;
		sun.getHorizontalRadiation(toa, direct, diffuse);
		if(param!=MeteoData::RSWR)
			value = (direct+diffuse)*solarIndex; //ISWR
		else
			value = (direct+diffuse)*albedo*solarIndex; //RSWR
	}

	return true; //all missing values could be filled
}

bool AllSkySWGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if(vecMeteo.empty()) return true;

	const double lat = vecMeteo.front().meta.position.getLat();
	const double lon = vecMeteo.front().meta.position.getLon();
	const double alt = vecMeteo.front().meta.position.getAltitude();
	if(lat==IOUtils::nodata || lon==IOUtils::nodata || alt==IOUtils::nodata) return false;
	sun.setLatLon(lat, lon, alt);

	bool all_filled = true;
	for(size_t ii=0; ii<vecMeteo.size(); ii++) {
		double &value = vecMeteo[ii](param);
		if(value == IOUtils::nodata) {
			const double ISWR=vecMeteo[ii](MeteoData::ISWR), RSWR=vecMeteo[ii](MeteoData::RSWR), HS=vecMeteo[ii](MeteoData::HS);
			double TA=vecMeteo[ii](MeteoData::TA), RH=vecMeteo[ii](MeteoData::RH), ILWR=vecMeteo[ii](MeteoData::ILWR);

			double albedo = .5;
			if(RSWR==IOUtils::nodata || ISWR==IOUtils::nodata) {
				if(HS!=IOUtils::nodata) //no big deal if we can not adapt the albedo
					albedo = (HS>=snow_thresh)? snow_albedo : soil_albedo;
			} else if(ISWR>0. && RSWR>0.) { //this could happen if the user calls this generator for a copy parameter, etc
				albedo = RSWR / ISWR;
				if(albedo>=1.) albedo=0.99;
				if(albedo<=0.) albedo=0.01;
			}

			if(TA==IOUtils::nodata || RH==IOUtils::nodata) {
				//set TA & RH so the reduced precipitable water will get an average value
				TA=274.98;
				RH=0.666;
				ILWR=IOUtils::nodata; //skip solarIndex correction
			}

			sun.setDate(vecMeteo[ii].date.getJulian(true), 0.);
			const double solarIndex = (ILWR!=IOUtils::nodata)? getSolarIndex(TA, RH, ILWR) : 1.;

			const double P=vecMeteo[ii](MeteoData::P);
			if(P==IOUtils::nodata)
				sun.calculateRadiation(TA, RH, albedo);
			else
				sun.calculateRadiation(TA, RH, P, albedo);

			double toa, direct, diffuse;
			sun.getHorizontalRadiation(toa, direct, diffuse);
			if(param!=MeteoData::RSWR)
				value = (direct+diffuse)*solarIndex; //ISWR
			else
				value = (direct+diffuse)*albedo*solarIndex; //RSWR
		}
	}

	return all_filled;
}

double AllSkySWGenerator::getSolarIndex(const double& ta, const double& rh, const double& ilwr)
{// this is based on Kartsen cloudiness, Dilley clear sky emissivity and Unsworth ILWR
//this means that this solar index is the ratio of iswr for clear sky on a horizontal
//surface and the measured iswr
	const double epsilon_clear = Atmosphere::Dilley_emissivity(rh, ta);
	const double ilwr_clear = Atmosphere::blkBody_Radiation(1., ta);

	double cloudiness = (ilwr/ilwr_clear - epsilon_clear) / (.84 * (1.-epsilon_clear));
	if(cloudiness>1.) cloudiness=1.;
	if(cloudiness<0.) cloudiness=0.;

	const double b1 = 0.75, b2 = 3.4;
	const double karsten_Si = 1. - (b1 * pow(cloudiness, b2));
	return karsten_Si;
}


void ESOLIPGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	if(vecArgs.size()>0) { //incorrect arguments, throw an exception
		throw InvalidArgumentException("Wrong number of arguments supplied for the "+algo+" generator", AT);
	}
}

bool ESOLIPGenerator::generate(const size_t& /*param*/, MeteoData& /*md*/)
{//HACK: modify prototype so we can get the full vector + the index of the replacement
	return false; //all missing values could be filled
}

//when we can not guarantee PSUM=0, we leave it at nodata. Therefore, it is highly recommended to
//run through a Cst=0 data generator afterward
bool ESOLIPGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (param!=MeteoData::PSUM)
		throw InvalidArgumentException("Trying to use "+algo+" generator on " + MeteoData::getParameterName(param) + " but it can only be applied to PSUM!!", AT);

	if (vecMeteo.empty()) return true;

	//Find first point that is not IOUtils::nodata
	size_t last_good = IOUtils::npos;
	for (size_t ii=0; ii<vecMeteo.size(); ii++){
		if (vecMeteo[ii](MeteoData::HS) != IOUtils::nodata){
			last_good = ii;
			break;
		}
	}

	if (last_good == IOUtils::npos) //can not find a good point to start
		return false;

	bool all_filled = (last_good>0)? false : true;
	for (size_t ii=last_good+1; ii<vecMeteo.size(); ii++) {
		const double HS_curr = vecMeteo[ii](MeteoData::HS);
		if(HS_curr==IOUtils::nodata) continue;

		const size_t start_idx = last_good+1;
		const double HS_prev = vecMeteo[last_good](MeteoData::HS);
		const double HS_delta = HS_curr - HS_prev;

		if (HS_delta>0.) {
			const double rho = newSnowDensity(vecMeteo[ii]);
			const double precip = HS_delta * rho; //in kg/m2 or mm
			ProcPSUMDistribute::SmartDistributePSUM(precip, start_idx, ii, param, vecMeteo);
		} else {
			all_filled = false;
		}

		last_good=ii;
	}

	return all_filled;
}

double ESOLIPGenerator::newSnowDensity(const MeteoData& md) const
{ //C. Zwart, "Significance of new-snow properties for snowcover development",
//master's thesis, 2007, Institute for Marine and Atmospheric Research, University of Utrecht, 78 pp.
	const double vw = max(2., md(MeteoData::VW));
	const double rh = md(MeteoData::RH);
	const double ta = md(MeteoData::TA) - Cst::t_water_triple_pt;
	const double beta01=3.28, beta1=0.03, beta02=-0.36, beta2=-0.75, beta3=0.3;

	double arg = beta01 + beta1*ta + beta2*asin(sqrt(rh)) + beta3*log10(vw);
	if (ta>=-14.)
		arg += beta02; // += beta2*ta;

	return min( max(30., pow(10., arg)), 250. ); //limit the density to the [30, 250] kg/m3 range
}


void PPhaseGenerator::parse_args(const std::vector<std::string>& vecArgs)
{
	const size_t nArgs = vecArgs.size();
	
	if (nArgs<1 || IOUtils::isNumeric(vecArgs[0]))
		throw InvalidArgumentException("Wrong arguments supplied to the "+algo+" generator. Please provide the method to use and its arguments!", AT);
	
	const std::string user_algo = IOUtils::strToUpper(vecArgs[0]);
	if (user_algo=="THRESH") {
		if (nArgs!=2)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" generator for the "+user_algo+" method", AT);
		IOUtils::convertString(fixed_thresh, vecArgs[1]);
		model = THRESH;
	} else if (user_algo=="RANGE") {
		if (nArgs!=3)
			throw InvalidArgumentException("Wrong number of arguments supplied to the "+algo+" generator for the "+user_algo+" method", AT);
		double range_thresh1, range_thresh2;
		IOUtils::convertString(range_thresh1, vecArgs[1]);
		IOUtils::convertString(range_thresh2, vecArgs[2]);
		if (range_thresh1==range_thresh2)
			throw InvalidArgumentException(algo+" generator, "+user_algo+" method: the two provided threshold must be different", AT);
		if (range_thresh1>range_thresh2) 
			std::swap(range_thresh1, range_thresh2);
		range_start = range_thresh1;
		range_norm = 1. / (range_thresh2-range_thresh1);
		model = RANGE;
	} else
		throw InvalidArgumentException("Unknown parametrization \""+user_algo+"\" supplied to the "+algo+" generator", AT);
}

bool PPhaseGenerator::generate(const size_t& param, MeteoData& md)
{
	double &value = md(param);
	if (value==IOUtils::nodata) {
		const double TA=md(MeteoData::TA);
		if (TA==IOUtils::nodata) return false;
		
		if (model==THRESH) {
			value = (TA>=fixed_thresh)? 1. : 0.;
		} else if (model==RANGE) {
			const double tmp_rainfraction = range_norm * (TA - range_start);
			value = (tmp_rainfraction>1)? 1. : (tmp_rainfraction<0.)? 0. : tmp_rainfraction;
		}
	}

	return true; //all missing values could be filled
}

bool PPhaseGenerator::generate(const size_t& param, std::vector<MeteoData>& vecMeteo)
{
	if (vecMeteo.empty()) return true;
	
	bool all_filled = true;
	for (size_t ii=0; ii<vecMeteo.size(); ii++) {
		if (!generate(param, vecMeteo[ii]))
			all_filled = false;
	}

	return all_filled;
}

} //namespace

