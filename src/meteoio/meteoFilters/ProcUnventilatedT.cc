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
#include <meteoio/meteoFilters/ProcUnventilatedT.h>
#include <cmath>

using namespace std;

namespace mio {

const double ProcUnventilatedT::dflt_albedo = .23;
const double ProcUnventilatedT::vw_thresh = 0.1; //wind speed threshold

ProcUnventilatedT::ProcUnventilatedT(const std::vector<std::string>& vec_args, const std::string& name)
                  : ProcessingBlock(name), usr_albedo(dflt_albedo),
                    usr_vw_thresh(IOUtils::nodata), nakamura(false)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcUnventilatedT::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	ovec = ivec;

	if (usr_vw_thresh!=IOUtils::nodata)
		filterTA(param, ovec);
	else
		correctTA(param, ovec);

}

void ProcUnventilatedT::filterTA(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	for (size_t ii=0; ii<ovec.size(); ii++) {
		const double vw = ovec[ii](MeteoData::VW);
		if (vw!=IOUtils::nodata && vw<usr_vw_thresh)
			ovec[ii](param) = IOUtils::nodata;
	}
}

void ProcUnventilatedT::correctTA(const unsigned int& param, std::vector<MeteoData>& ovec) const
{
	if (nakamura) {
		for (size_t ii=0; ii<ovec.size(); ii++) {
			double& ta = ovec[ii](param);
			if (ta == IOUtils::nodata) continue; //preserve nodata values

			double iswr = ovec[ii](MeteoData::ISWR);
			const double rswr = ovec[ii](MeteoData::RSWR);
			double vw = ovec[ii](MeteoData::VW);
			const double hs = ovec[ii](MeteoData::HS);

			if (hs!=IOUtils::nodata && iswr==IOUtils::nodata && rswr!=IOUtils::nodata) {
				if (hs>snow_thresh) iswr = rswr / snow_albedo;
				else iswr = rswr / soil_albedo;
			}

			if (iswr==IOUtils::nodata || vw==IOUtils::nodata)
				continue;

			if (vw<vw_thresh) vw = vw_thresh; //this should be around the minimum measurable wind speed on regular instruments
			const double rho = 1.2; // in kg/m3
			const double Cp = 1004.;
			const double X = iswr / (rho*Cp*ta*vw);
			if (X<1e-4) continue; //the correction does not work well for small X values

			const double C0 = 0.13;
			const double C1 = 373.40;
			const double RE = C0 + C1*X;
			ta -= RE; //substracting the radiative error
		}
	} else { //Huwald
		for (size_t ii=0; ii<ovec.size(); ii++) {
			double& ta = ovec[ii](param);
			if (ta == IOUtils::nodata) continue; //preserve nodata values

			const double iswr = ovec[ii](MeteoData::ISWR);
			double rswr = ovec[ii](MeteoData::RSWR);
			double vw = ovec[ii](MeteoData::VW);
			const double hs = ovec[ii](MeteoData::HS);

			if (hs!=IOUtils::nodata && rswr==IOUtils::nodata && iswr!=IOUtils::nodata) {
				if (hs>snow_thresh) rswr = iswr * snow_albedo;
				else rswr = iswr * soil_albedo;
			}

			if (rswr==IOUtils::nodata || vw==IOUtils::nodata)
				continue;

			if (vw<vw_thresh) vw = vw_thresh; //this should be around the minimum measurable wind speed on regular instruments
			const double rho = 1.2; // in kg/m3
			const double Cp = 1004.;
			const double X = rswr / (rho*Cp*ta*vw) * 1e3;
			if (X<0.05) continue; //the correction does not work well for small X values

			const double RE = 3.1 * sqrt(X);
			ta -= RE; //substracting the radiative error
		}
	}
}

void ProcUnventilatedT::parse_args(std::vector<std::string> vec_args)
{
	const size_t nrArgs = vec_args.size();

	//no args -> simple correction with default albedo
	if (nrArgs==0) return;

	const bool arg0_is_num = IOUtils::isNumeric(vec_args[0]);

	if (nrArgs==1) { //either nakamura or default albedo
		if ( arg0_is_num ) {
			IOUtils::convertString(usr_albedo, vec_args[0]);
		} else {
			const string strArg( IOUtils::strToUpper(vec_args[0]) );
			if (strArg=="NAKAMURA")
				nakamura = true;
			else if (strArg=="HUWALD")
				nakamura = false;
			else if (strArg=="SUPPR")
				throw InvalidArgumentException("Invalid use of the SUPPR option for filter " + getName(), AT);
			else
				throw InvalidArgumentException("Unknown argument \""+strArg+"\" for filter " + getName(), AT);
		}
	} else if (nrArgs==2) {
		const string strArg = (!arg0_is_num)? IOUtils::strToUpper(vec_args[0]) : IOUtils::strToUpper(vec_args[1]);
		double dblArg;
		if (arg0_is_num)
			IOUtils::convertString(dblArg, vec_args[0]);
		else
			IOUtils::convertString(dblArg, vec_args[1]);

		if (strArg=="NAKAMURA") {
			nakamura = true;
			usr_albedo = dblArg;
		} else if (strArg=="HUWALD") {
			nakamura = false;
			usr_albedo = dblArg;
		} else if (strArg=="SUPPR") {
			usr_vw_thresh = dblArg;
		} else
			throw InvalidArgumentException("Unknown argument \""+strArg+"\" for filter " + getName(), AT);
	} else
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);

	if (usr_albedo<0. || usr_albedo>1.)
		throw InvalidArgumentException("Albedo value should be between 0 and 1 for filter " + getName(), AT);
}

} //end namespace
