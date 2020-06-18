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
 * @file ReSolver1d.h
 * @version 10.02
 */

#ifndef RESOLVER1D_H
#define RESOLVER1D_H

#include "../DataClasses.h"
#include "../SnowpackConfig.h"
#include <meteoio/MeteoIO.h>

#include <string.h>
/**
 * @class ReSolver1d
 * @version 10.02
 * @author Nander Wever
 * @bug Prone to bugs at any changes! Be aware!
 * @brief This module contains the solver for the 1d Richards Equation for the 1d snowpack model
 */
class ReSolver1d {

	public:
		ReSolver1d(const SnowpackConfig& cfg);	// Class constructor
		void SolveRichardsEquation(SnowStation& Xdata, SurfaceFluxes& Sdata);

		double surfacefluxrate;		// Surfacefluxrate for solving RE. It is either surface of snow, in case of snowpack and solving RE for snow, or surface of soil, when no snowpack and/or solving RE only for soil.
		double soilsurfacesourceflux;	// Soilsurfacesourceflux for solving RE. This is used when we use RE for snow AND there is a snowpack AND the lowest snow element is removed.


	private:
		std::string variant;

		//To prevent string comparisons, we define an enumerated list:
		enum watertransportmodels{UNDEFINED, BUCKET, NIED, RICHARDSEQUATION};
		//Soil types
		enum SoilTypes{ORGANIC, CLAY, CLAYLOAM, LOAM, LOAMYSAND, SAND, SANDYCLAY, SANDYCLAYLOAM, SANDYLOAM, SILT, SILTYCLAY, SILTYCLAYLOAM, SILTLOAM, WFJGRAVELSAND};
		//Hydraulic conductivity parameterizations
		enum K_Parameterizations{SHIMIZU, CALONNE};
		//K_Average types
		enum K_AverageTypes{ARITHMETICMEAN, GEOMETRICMEAN, HARMONICMEAN, MINIMUMVALUE, UPSTREAM};
		//Van genuchten model types
		enum VanGenuchten_ModelTypesSnow{YAMAGUCHI2012, YAMAGUCHI2010, YAMAGUCHI2010_ADAPTED, DAANEN};
		//Solvers
		enum SOLVERS{DGESVD, DGTSV, TDMA};
		//Boundary conditions
		enum BoundaryConditions{DIRICHLET, NEUMANN, LIMITEDFLUXEVAPORATION, LIMITEDFLUXINFILTRATION, LIMITEDFLUX, WATERTABLE, FREEDRAINAGE, GRAVITATIONALDRAINAGE, SEEPAGEBOUNDARY};
		
		
		watertransportmodels iwatertransportmodel_snow, iwatertransportmodel_soil;

		std::string watertransportmodel_snow;
		std::string watertransportmodel_soil;
		BoundaryConditions BottomBC;				//Bottom boundary condition (recommended choice either DIRICHLET with saturation (lower boundary in water table) or FREEDRAINAGE (lower boundary not in water table))
		K_AverageTypes K_AverageType;				//Implemented choices: ARITHMETICMEAN (recommended), HARMONICMEAN, GEOMETRICMEAN, MINIMUMVALUE, UPSTREAM

		double sn_dt;
		bool useSoilLayers, water_layer;


		// Van Genuchten functions
		double fromTHETAtoH(double theta, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e, double h_d);
		double fromTHETAtoHforICE(double theta, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e, double h_d, double theta_i);
		double fromHtoTHETA(double h, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e);
		double fromHtoTHETAforICE(double h, double theta_r, double theta_s, double alpha, double m, double n, double Sc, double h_e, double theta_i);
		double AirEntryPressureHead(double MaximumPoreSize, double Temperature);
		void SetSoil(SoilTypes type, double *theta_r, double *theta_s, double *alpha, double *m, double *n, double *ksat, double *he);

		// Solvers
		int TDMASolver (int n, double *a, double *b, double *c, double *v, double *x);
		int pinv(int m, int n, int lda, double *a);
};
#endif //End of WaterTransport.h
