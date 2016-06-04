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
#include <meteoio/meteoFilters/ProcUndercatch_WMO.h>
#include <meteoio/meteoLaws/Atmosphere.h>
#include <cmath>

using namespace std;

namespace mio {

const double ProcUndercatch_WMO::Tsnow_WMO=-2., ProcUndercatch_WMO::Train_WMO=2.; //WMO values from Yan et al (2001)

ProcUndercatch_WMO::ProcUndercatch_WMO(const std::vector<std::string>& vec_args, const std::string& name)
                   : ProcessingBlock(name), type(cst),
                     factor_snow(1.3), factor_mixed(1.1), Tsnow(Tsnow_WMO), Train(Train_WMO)
{
	parse_args(vec_args);
	properties.stage = ProcessingProperties::first; //for the rest: default values
}

void ProcUndercatch_WMO::process(const unsigned int& param, const std::vector<MeteoData>& ivec,
                        std::vector<MeteoData>& ovec)
{
	if(param!=MeteoData::PSUM)
		throw InvalidArgumentException("Trying to use "+getName()+" filter on " + MeteoData::getParameterName(param) + " but it can only be applied to precipitation!!" + getName(), AT);
	ovec = ivec;

	for (size_t ii=0; ii<ovec.size(); ii++){
		double& tmp = ovec[ii](param);
		double VW = ovec[ii](MeteoData::VW);
		if(VW!=IOUtils::nodata) VW = Atmosphere::windLogProfile(VW, 10., 2.); //impact seems minimal
		double t = ovec[ii](MeteoData::TA);
		if(t==IOUtils::nodata) continue; //we MUST have air temperature in order to filter
		t=IOUtils::K_TO_C(t); //t in celsius
		precip_type precip = (t<=Tsnow)? snow : (t>=Train)? rain : mixed;

		//We don't use Tmax, Tmin, Tmean but only the current temperature instead
		if (tmp == IOUtils::nodata || tmp==0.) {
			continue; //preserve nodata values and no precip
		} else if(type==cst) {
			if(precip==snow) tmp *= factor_snow;
			if(precip==mixed) tmp *= factor_mixed;
		} else if(type==nipher) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(precip==snow) k=100.-0.44*VW*VW-1.98*VW;
			if(precip==mixed) {
				k=97.29-3.18*VW+0.58*t-0.67*t; //Tmax, Tmin
			}
			tmp *= 100./k;
		} else if(type==tretyakov) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(VW>8.5) VW=8.5; //the fits have been calibrated until 8.5 m/s
			if(precip==snow) k=103.11-8.67*VW+0.30*t; //Tmax
			if(precip==mixed) k=96.99-4.46*VW+0.88*t+0.22*t; //Tmax, Tmin
			tmp *= 100./k;
		} else if(type==us8sh) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(precip==snow) k=exp(4.61-0.04*pow(VW, 1.75));
			if(precip==mixed) k=101.04-5.62*VW;
			tmp *= 100./k;
		} else if(type==us8unsh) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(precip==snow) k=exp(4.61-0.16*pow(VW, 1.28));
			if(precip==mixed) k=100.77-8.34*VW;
			tmp *= 100./k;
		} else if(type==rt3_jp) {
			if(VW==IOUtils::nodata) continue;
			const double rh = ovec[ii](MeteoData::RH);
			const double alt = ovec[ii].meta.position.getAltitude();
			double k=100.;
			if(rh!=IOUtils::nodata && alt!=IOUtils::nodata) {
				const double t_wb = IOUtils::K_TO_C(Atmosphere::wetBulbTemperature(ovec[ii](MeteoData::TA), rh, alt));
				double ts_rate;
				if(t_wb<1.1) ts_rate = 1. - .5*exp(-2.2*pow(1.1-t_wb, 1.3));
				else ts_rate = .5*exp(-2.2*pow(t_wb-1.1, 1.3));
				if(ts_rate>.5) precip=snow; else precip=mixed;
			}
			if(precip==snow) k=100. / (1.+.346*VW);
			if(precip==mixed) k=100. / (1.+.0856*VW);
			tmp *= 100./k;
		} else if(type==cspg) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(precip==snow) k=100.*exp(-0.056*VW); //VW in 0 - 6.2
			if(precip==rain) k=100.*exp(-0.041*VW); //VW in 0 - 7.3
			if(precip==mixed) {
				const double Ksnow = 100.*exp(-0.056*VW);
				const double Krain = 100.*exp(-0.041*VW);
				const double td = (t<-2.)? -2. : (t>2.)? 2. : t;
				k = Ksnow - (Ksnow-Krain)*(td+2.)/4.; //Tmean
			}
			tmp *= 100./k;
		} else if(type==geonorsh) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(precip==snow) k=100.*exp(-0.135*VW); //VW in 0 - 6
			if(precip==rain) k=100.*exp(-0.113*VW); //VW in 0 - 5
			if(precip==mixed) {
				const double Ksnow = 100.*exp(-0.135*VW);
				const double Krain = 100.*exp(-0.113*VW);
				const double td = (t<-2.)? -2. : (t>2.)? 2. : t;
				k = Ksnow - (Ksnow-Krain)*(td+2.)/4.; //Tmean
			}
			tmp *= 100./k;
		} else if(type==hellmann) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(precip==snow) k=100.+1.13*VW*VW-19.45*VW;
			if(precip==mixed) k=96.63+0.41*VW*VW-9.84*VW+5.95*t; //Tmean
			tmp *= 100./k;
		}  else if(type==hellmannsh) {
			if(VW==IOUtils::nodata) continue;
			double k=100.;
			if(precip==snow) k=100.+0.72*VW*VW-13.74*VW;
			if(precip==mixed) k=101.319+0.524*VW*VW-6.42*VW;
			tmp *= 100./k;
		}
	}
}

void ProcUndercatch_WMO::parse_args(std::vector<std::string> filter_args)
{
	if (filter_args.empty())
		throw InvalidArgumentException("Wrong number of arguments for filter " + getName(), AT);

	for(size_t ii=0; ii<filter_args.size(); ii++) {
		IOUtils::toLower(filter_args[ii]);
	}

	if(filter_args[0]=="cst") {
		type=cst;
		if(filter_args.size() < 3 || filter_args.size() > 5 || filter_args.size() == 4)
			throw InvalidArgumentException("Wrong number of arguments for filter "+getName()+" with rain gauge type \"cst\"", AT);
		IOUtils::convertString(factor_snow, filter_args[1]);
		IOUtils::convertString(factor_mixed, filter_args[2]);
		if(filter_args.size()==5) {
			IOUtils::convertString(Tsnow, filter_args[3]);
			Tsnow = IOUtils::K_TO_C(Tsnow);
			IOUtils::convertString(Train, filter_args[4]);
			Train = IOUtils::K_TO_C(Train);
		}
	} else if(filter_args[0]=="nipher") {
		type=nipher;
	} else if(filter_args[0]=="tretyakov") {
		type=tretyakov;
	} else if(filter_args[0]=="us8sh") {
		type=us8sh;
	} else if(filter_args[0]=="us8unsh") {
		type=us8unsh;
	} else if(filter_args[0]=="rt3_jp") {
		type=rt3_jp;
	} else if(filter_args[0]=="cspg") {
		type=cspg;
	} else if(filter_args[0]=="geonorsh") {
		type=geonorsh;
	} else if(filter_args[0]=="hellmann") {
		type=hellmann;
	} else if(filter_args[0]=="hellmannsh") {
		type=hellmannsh;
	} else {
		throw InvalidArgumentException("Rain gauge type \""+ filter_args[0] +"\" unknown for filter "+getName(), AT);
	}
}

} //end namespace
