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
 * @file WaterTransport.h
 * @version 10.02
 */

#ifndef WATERTRANSPORT_H
#define WATERTRANSPORT_H

#include <snowpack/DataClasses.h>
#include <snowpack/snowpackCore/ReSolver1d.h>

#include <meteoio/MeteoIO.h>

/**
 * @class WaterTransport
 * @version 10.02
 * @bug Prone to bugs at any changes! Be aware!
 * @brief This module contains water transport routines for the 1d snowpack model
 */
class WaterTransport {

	public:
		WaterTransport(const SnowpackConfig& cfg);
		void compTransportMass(const CurrentMeteo& Mdata, const double& ql, SnowStation& Xdata, SurfaceFluxes& Sdata);

	private:
		//The following 3 functions are used in WaterTransport model "NIED"
		double BisFunc(const double X, const double P[]);
		double Bisection(const double minval, const double maxval, double P[]);
		void KHCalcNaga(const double RG, const double Dens, double ThR, const double WatCnt, const double SatuK, double &Rh, double &Rk);

		void compSurfaceSublimation(const CurrentMeteo& Mdata, double ql, SnowStation& Xdata, SurfaceFluxes& Sdata);

		void mergingElements(SnowStation& Xdata, SurfaceFluxes& Sdata);

		void adjustDensity(SnowStation& Xdata);

		void transportWater(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata);

		ReSolver1d RichardsEquationSolver1d;

		std::string variant;

		//To prevent string comparisons, we define an enumerated list:
		enum watertransportmodels{UNDEFINED, BUCKET, NIED, RICHARDSEQUATION};
		watertransportmodels iwatertransportmodel_snow, iwatertransportmodel_soil;

		std::string watertransportmodel_snow;
		std::string watertransportmodel_soil;
		double sn_dt;
		double hoar_thresh_rh, hoar_thresh_vw, hoar_thresh_ta;
		double hoar_density_buried, hoar_density_surf, hoar_min_size_buried;
		double minimum_l_element;
		bool useSoilLayers, water_layer, jam;
};
#endif //End of WaterTransport.h

