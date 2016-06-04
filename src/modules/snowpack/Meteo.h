/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of Snowpack.
    Snowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snowpack.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @file Meteo.h
 * @author Charles Fierz
 * @version 9.x
 * @date -
 */

#ifndef __METEO_H__
#define __METEO_H__

#include <meteoio/MeteoIO.h>

#include <snowpack/SnowpackConfig.h>
#include <snowpack/Constants.h>
#include <snowpack/Hazard.h>
#include <snowpack/Utils.h>
#include <snowpack/Laws_sn.h>
#include <snowpack/snowpackCore/Snowpack.h>
#include <snowpack/snowpackCore/Canopy.h>

class Meteo {

	public:
		typedef enum {
			RICHARDSON,  ///< Simplified Richardson number stability correction
			MONIN_OBUKHOV, ///< Standard MO iteration with Paulson and Stearns C. and Weidner G., <i>"sensible and latent heat flux estimates in antarctica"</i>, Antarctic meteorology and climatology: studies based on automatic weather stations, Antarctic Research Series, <b>61</b>, pp 190--138, 1993
			NEUTRAL_MO  ///< Assume neutral MO stratification
		} ATM_STABILITY;

		Meteo(const SnowpackConfig& i_cfg);

		static void projectPrecipitations(const double& SlopeAngle, double& precips, double& hs);
		static bool compHSrate(CurrentMeteo& Mdata, const SnowStation& vecXdata, const double& hs_a3hl6);
		void compMeteo(CurrentMeteo &Mdata, SnowStation &Xdata, const bool& runCanopyModel);
		static void compRadiation(const SnowStation &station, mio::SunObject &sun, SnowpackConfig &cfg, CurrentMeteo &Mdata);
		static void radiationOnSlope(const SnowStation &sector, const mio::SunObject &sun, CurrentMeteo &Mdata, SurfaceFluxes &surfFluxes);
		void setStability(const ATM_STABILITY& i_stability);
		ATM_STABILITY getStability() const;

 	private:
		void MicroMet(const SnowStation& Xdata, CurrentMeteo& Mdata, const bool& adjust_VW_height=true) const;
		static double getParameterAverage(mio::IOManager& io, const mio::MeteoData::Parameters& param,
		                                  const mio::Date& current_date, const int& time_span, const int& increment);
		static void RichardsonStability(const double& ta_v, const double& t_surf_v, const double& zref,
		                                const double& vw, const double& z_ratio, double &ustar, double &psi_s);
		static void MOStability(const double& ta_v, const double& t_surf_v, const double& t_surf,
		                        const double& zref, const double& vw, const double& z_ratio, double &ustar,
		                        double &psi_s, double &psi_m);

		Canopy canopy;
		double roughness_length, height_of_wind_value;
		ATM_STABILITY stability;
		bool research_mode, useCanopyModel, alpine3d;
};

#endif //END of Meteo.h
