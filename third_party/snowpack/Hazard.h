/* **********************************************************************************************/
/*                                        stand-alone                                          */
/*                               Derived from RESEARCH VERSION 9.0                             */
/* **********************************************************************************************/
/* **********************************************************************************/
/*  Copyright WSL Institute for Snow and Avalanche Research    SLF-DAVOS           */
/* **********************************************************************************/
/* This file is part of Snowpack.
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
 * @file Hazard.h
 * @version 10.02
 * This module contains the hazard computation routines and structures
*/

#ifndef HAZARD_H
#define HAZARD_H

#include "DataClasses.h"
#include <meteoio/MeteoIO.h>
#include <vector>

/** @brief 
 *
 * @ingroup postprocessing
 */
class Hazard {
	public:
		Hazard(const SnowpackConfig& cfg, const double duration);

		void initializeHazard(std::vector<double>& vecDrift, double slope_angle,
		                      std::vector<ProcessDat>& Hdata, std::vector<ProcessInd>& Hdata_ind);

		void getHazardDataMainStation(ProcessDat& Hdata, ProcessInd& Hdata_ind,
                                      ZwischenData& Zdata, const double& newDrift, const bool stationDriftIndex,
                                      const SnowStation& Xdata, const CurrentMeteo& Mdata, const SurfaceFluxes& Sdata);

		void getHazardDataSlope(ProcessDat& Hdata, ProcessInd& Hdata_ind,
		                        std::vector<double>& drift24, const double& newDrift, const SnowStation& Xdata,
		                        const bool luvDriftIndex, const bool north, const bool south);

		static const double typical_slope_length, wind_slab_density;
		static const double minimum_drift, maximum_drift;

	private:
		enum ActVec {
			noAction=0,     // neither shift nor overwite index values
			overwrite,      // overwrite <index>[0] w/o shifting
			pushOverwrite   // push vector values and overwrite <index>[0]
		};

		void actOnVector(std::vector<double>& oldVector, const double& newValue, const ActVec& action);

		double compDriftIndex(std::vector<double>& vecDrift, const double& drift, const double& rho,
		                      const unsigned int& nHours, const double& slope_angle, const ActVec& action);

		void getDriftIndex(ProcessDat& Hdata, ProcessInd& Hdata_ind,
		                   std::vector<double>& vecDrift, const double& newDriftValue, const double slope_angle);

		double compHoarIndex(std::vector<double> &oldHoar, const double& newHoar,
                             const unsigned int& nHours, const ActVec& action);

		double compDewPointDeficit(double TA, double TSS, double RH);

		void compMeltFreezeCrust(const SnowStation& Xdata, ProcessDat& Hdata, ProcessInd& Hdata_ind);

		bool research_mode, enforce_measured_snow_heights, force_rh_water;
		unsigned int nHz, hazard_steps_between;
		double sn_dt;
		double hoar_density_surf, hoar_min_size_surf;
};

#endif
