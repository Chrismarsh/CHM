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
#ifndef __SALTATION_H__
#define __SALTATION_H__

#include <snowpack/Utils.h>
#include <meteoio/MeteoIO.h>
#include <cmath>

/**
 * @brief class Saltation
 * @author Michael Lehning \n Judith Doorschot
 * @version 9.x
 * @date    -
 * @bug     -
 * @brief This module contains the saltation model of Judith.
 * - 24.08.2000: The new very complicated model. Judith says that it is my fault, if it is wrong. I hope not ..... \n
 *               Make a separate routine to model saltation only.
 * - 26.04.2001: Finally, on the Friday evening before the Swiss Bike Masters event, where Michael was
 *               supposed to start, he also started to implement the last version of Judith's saltation
 *               model. The GRID man Tuan Anh (Nguyen) had arrived and smiled. \n
 * The Gaudergrat experiment GAUDEX was in good shape and almost everything was up and
 * running.
 */
class Saltation {
	public:
		Saltation(const SnowpackConfig& i_cfg);

		bool compSaltation(const double& tauS, const double& tau_th, const double& slope_angle, const double& dg,
		                   double& massflux, double& c_salt);

		static const double karman;
		static const double z0_salt;

	private:
		double sa_vw(const double& z, const double& tauA, const double& tauS, const double& z0,
                     const double& u_start, const double& slope_angle);
		double sa_vw2(const double& z, const double& tauA, const double& tauS, const double& z0,
                      const double& u_start, const double& slope_angle);

		bool sa_Traject(const double& u0, const double& angle_e_rad, const double& slope_angle, const double& dg,
		                const double& tauA, const double& tauS, const double& z0,
		                double& ubar, double& u_i, double& angle_i_rad, double& t_i, double& z_max);

		double sa_MassFlux(const double& z0, const double& tauS, const double& tauA, const double& slope_angle,
		                   const double& dg, const double& tau_th, double& z_max, double& ubar, double& cs);

		double sa_AeroEntrain(const double& z0, const double& tauS, const double& slope_angle, const double& dg,
		                      const double& tau_th, double& flux, double& z_max, double& ubar, double& cs);

		int sa_TestSaltation(const double& z0, const double& tauS, const double& tauA, const double& slope_angle,
		                     const double& dg, const double& tau_th, double& z_max, double& ubar);

		std::string saltation_model;
		static const double hs_frac, elas, angle_ej, ratio_ve_ustar, salt_height;
		static const int strong, weak;
};

#include <snowpack/Constants.h>

#endif

