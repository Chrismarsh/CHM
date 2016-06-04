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

#ifndef __SNOWDRIFT_H__
#define __SNOWDRIFT_H__

#include <meteoio/MeteoIO.h>
#include <snowpack/Constants.h>
#include <cmath>
#include <snowpack/Saltation.h>
#include <snowpack/snowpackCore/Snowpack.h>
#include <vector>

class Saltation;

/**
 * @class SnowDrift
 * @brief This class contains the computation of local snow drift and the associated erosion
 * @version 10.02
 */
class SnowDrift {

	public:
		SnowDrift(const SnowpackConfig& i_cfg);

		void compSnowDrift(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& cumu_psum);

		static const double schmidt_drift_fudge;

 	private:
		double compMassFlux(const ElementData& Edata, const double& ustar, const double& slope_angle);

		Saltation saltation; // The saltation model used
		bool enforce_measured_snow_heights, snow_redistribution, snow_erosion; // Will be read from cfg object
		bool alpine3d; ///< triggers various tricks for Alpine3D (including reducing the number of warnings)
		double sn_dt;        //Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
		int nSlopes;
		static const bool msg_erosion;
}; //End class SnowDrift

#endif //#ifndef __SNOWDRIFT_H__

