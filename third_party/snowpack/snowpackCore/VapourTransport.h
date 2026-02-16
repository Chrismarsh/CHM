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
 * @file VapourTransport.h
 */

#ifndef VAPOURTRANSPORT_H
#define VAPOURTRANSPORT_H

#include "../Constants.h"
#include "../DataClasses.h"
#include "../Laws_sn.h"
#include "ReSolver1d.h"
#include "WaterTransport.h"
#include "Snowpack.h"
#include "PhaseChange.h"
#include "../Meteo.h"
#include "../Utils.h"
#include "Solver.h"
#include "../Constants.h"
#include "../Laws_sn.h"
#include "../SnowDrift.h"
#include "Metamorphism.h"

#include <meteoio/MeteoIO.h>


/**
 * @class VapourTransport
 * @version 1.0
 * @brief This module contains water vapour transport routines for the 1d snowpack model
 */
class VapourTransport : public WaterTransport {
	public:
		VapourTransport(const SnowpackConfig& cfg);
		void compTransportMass(const CurrentMeteo& Mdata, double& ql, SnowStation& Xdata, SurfaceFluxes& Sdata);

	private:
		bool compDensityProfile(const CurrentMeteo& Mdata, SnowStation& Xdata,
                                std::vector<double>& hm_, std::vector<double>& as_,
                                const std::vector<double>& D_el, std::vector<double>& oldVaporDenNode);
		void compSurfaceSublimation(const CurrentMeteo& Mdata, double& ql, SnowStation& Xdata, SurfaceFluxes& Sdata);
		void LayerToLayer(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& ql);
		double dRhov_dT(const double Tem);

		ReSolver1d RichardsEquationSolver1d;

		std::string variant;

		watertransportmodels iwatertransportmodel_snow, iwatertransportmodel_soil;

		std::string watertransportmodel_snow;
		std::string watertransportmodel_soil;
		double sn_dt, timeStep, waterVaporTransport_timeStep;
		const static double VapourTransport_timeStep;
		double hoar_thresh_rh, hoar_thresh_vw, hoar_thresh_ta;
		bool useSoilLayers, water_layer;

		bool enable_vapour_transport;
		double diffusionScalingFactor_, height_of_meteo_values;
		bool adjust_height_of_meteo_values;

		const static double f;

		bool waterVaporTransport_timeStepAdjust;
};
#endif // End of VapourTransport.h}
