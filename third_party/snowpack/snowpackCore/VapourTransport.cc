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


#include <snowpack/snowpackCore/VapourTransport.h>
#include <snowpack/vanGenuchten.h>
#include <snowpack/snowpackCore/Snowpack.h>
#include <snowpack/Constants.h>

// MeteoIO constants
#include <meteoio/meteoLaws/Meteoconst.h>

#include <assert.h>
#include <sstream>
#include <errno.h>

//Eigen, note we temporarily disable Effective C++ warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#include <meteoio/thirdParty/Eigen/Dense>
#include <meteoio/thirdParty/Eigen/Sparse>
#include <meteoio/thirdParty/Eigen/IterativeLinearSolvers>
#include <meteoio/thirdParty/Eigen/SparseQR>
#include <meteoio/thirdParty/Eigen/SparseCholesky>
#include <meteoio/thirdParty/Eigen/SparseLU>
#include <meteoio/thirdParty/Eigen/Core>

typedef Eigen::Triplet<double> Trip;
#pragma GCC diagnostic pop

using namespace mio;
using namespace std;
using namespace Eigen;


/**
 * @page water_vapor_transport Water Vapor Transport
 *
 */

//vapour_transport_implicit_factor: 1 is fully implicit, 0 is fully explicit, 0.5 is Crank-Nicolson
const double VapourTransport::f = 1.;
const double VapourTransport::VapourTransport_timeStep = 60.;	// Only used when f < 1 !!

VapourTransport::VapourTransport(const SnowpackConfig& cfg)
               : WaterTransport(cfg), RichardsEquationSolver1d(cfg, false), variant(),
                 iwatertransportmodel_snow(BUCKET), iwatertransportmodel_soil(BUCKET), watertransportmodel_snow("BUCKET"), watertransportmodel_soil("BUCKET"),
                 sn_dt(IOUtils::nodata), timeStep(IOUtils::nodata), waterVaporTransport_timeStep(IOUtils::nodata),
                 hoar_thresh_rh(IOUtils::nodata), hoar_thresh_vw(IOUtils::nodata), hoar_thresh_ta(IOUtils::nodata),
                 useSoilLayers(false), water_layer(false), enable_vapour_transport(false),
                 diffusionScalingFactor_(1.0), height_of_meteo_values(0.), adjust_height_of_meteo_values(true), waterVaporTransport_timeStepAdjust(false)
{
	cfg.getValue("VARIANT", "SnowpackAdvanced", variant);

	// Defines whether soil layers are used
	cfg.getValue("SNP_SOIL", "Snowpack", useSoilLayers);

	//To build a thin top rain-water layer over a thin top ice layer, rocks, roads etc.
	cfg.getValue("WATER_LAYER", "SnowpackAdvanced", water_layer);

	/**
	 * @brief No surface hoar will form for rH above threshold (1)
	 * - Original calibration with the 98/99 data set: 0.9
	 * - r141: HOAR_THRESH_RH set to 0.9
	 * - r719: HOAR_THRESH_RH set to 0.97
	 */
	cfg.getValue("HOAR_THRESH_RH", "SnowpackAdvanced", hoar_thresh_rh);

	/**
	 * @brief No surface hoar will form at wind speeds above threshold (m s-1)
	 * - Original calibration with the 98/99 data set: 3.5
	 * - r141: HOAR_THRESH_VW set to 3.0
	 * - r242: HOAR_THRESH_VW set to 3.5
	 */
	cfg.getValue("HOAR_THRESH_VW", "SnowpackAdvanced", hoar_thresh_vw);

	/**
	 * @brief No surface hoar will form at air temperatures above threshold (m s-1)
	 * - Originaly, using THRESH_RAIN
	 * - r787: HOAR_THRESH_TA set to 1.2
	 */
	cfg.getValue("HOAR_THRESH_TA", "SnowpackAdvanced", hoar_thresh_ta);

	//Calculation time step in seconds as derived from CALCULATION_STEP_LENGTH
	const double calculation_step_length = cfg.get("CALCULATION_STEP_LENGTH", "Snowpack");
	sn_dt = M_TO_S(calculation_step_length);

	//Vapour transport settings
	cfg.getValue("ENABLE_VAPOUR_TRANSPORT", "SnowpackAdvanced", enable_vapour_transport);
	if (enable_vapour_transport) {
		// the water vapor subtime step
		// If not using fully implicit scheme
		if (f < 1.0) {
			waterVaporTransport_timeStepAdjust = true;
			waterVaporTransport_timeStep = std::min(sn_dt, VapourTransport_timeStep);
		} else {
			// Implicit: time step defaults to SNOWPACK time step
			waterVaporTransport_timeStep = sn_dt;
		}
	}

	// Water transport model snow
	cfg.getValue("WATERTRANSPORTMODEL_SNOW", "SnowpackAdvanced", watertransportmodel_snow);
	iwatertransportmodel_snow = UNDEFINED;
	if (watertransportmodel_snow == "BUCKET") {
		iwatertransportmodel_snow = BUCKET;
	} else if (watertransportmodel_snow == "NIED") {
		iwatertransportmodel_snow = NIED;
	} else if (watertransportmodel_snow == "RICHARDSEQUATION") {
		iwatertransportmodel_snow = RICHARDSEQUATION;
	}

	// Water transport model soil
	cfg.getValue("WATERTRANSPORTMODEL_SOIL", "SnowpackAdvanced", watertransportmodel_soil);
	iwatertransportmodel_soil = UNDEFINED;
	if (watertransportmodel_soil == "BUCKET") {
		iwatertransportmodel_soil = BUCKET;
	} else if (watertransportmodel_soil == "NIED") {
		iwatertransportmodel_soil = NIED;
	} else if (watertransportmodel_soil == "RICHARDSEQUATION") {
		iwatertransportmodel_soil = RICHARDSEQUATION;
	}

	cfg.getValue("HEIGHT_OF_METEO_VALUES", "Snowpack", height_of_meteo_values);
	cfg.getValue("ADJUST_HEIGHT_OF_METEO_VALUES", "SnowpackAdvanced", adjust_height_of_meteo_values);
}

/**
 * @brief The mass transport procedure, which serves as the primary function, is invoked from Snowpack::runSnowpackModel. \n
 * NOTES:
 * -#   It is worth noting that the solver is highly stable with default parameters. Specifically, VAPOUR_TRANSPORT_TIMESTEP is set to 60 seconds,
 *      while the SNOWPACK simulation time step is 15 minutes, and VAPOUR_TRANSPORT_IMPLICIT_FACTOR is set to 1. The latter factor determines whether
 *      the equation is discretized in full implicit, full explicit, or a combination of the two, with a value of 1 indicating full implicit,
 *      and a value of 0.5 indicating Crank-Nicolson. In the case of convergence issues, reducing the height of the new-snow element controlled
 *      by HEIGHT_NEW_ELEM (in the .ini config file) is recommended. For sea-ice simulations, choosing BUCKET for the water transport scheme is advised
 *      if convergence issues arise.  \n
 * -#   If there is no soil or snow present, vapor transport will be bypassed.
 * -#   If vapor transport enabled, ql is only used in vaportransport for mass tranport on top. See WaterTransport::compTransportMass. \n
 * @author Mahdi Jafari
 * @param Xdata
 * @param ql Latent heat flux (W m-2)
 * @param Sdata
 * @param Mdata
 */
void VapourTransport::compTransportMass(const CurrentMeteo& Mdata, double& ql,
                                       SnowStation& Xdata, SurfaceFluxes& Sdata)
{

	// First, consider no soil with no snow on the ground
	if (!useSoilLayers && Xdata.getNumberOfNodes() == Xdata.SoilNode+1) {
		return;
	}

	try {
		LayerToLayer(Mdata, Xdata, Sdata, ql);
		WaterTransport::adjustDensity(Xdata);
	} catch(const exception&) {
		prn_msg( __FILE__, __LINE__, "err", Mdata.date, "Error in transportVapourMass()");
		throw;
	}
}


/**
 * @brief This function prepares everything to solve the transient-diffusive vapor transport with phase change:  \n
 * NOTES:
 * -#   Initially, the model employs the complete latent heat flux (ql) to alter the mass of the uppermost components.
 *      It is crucial to note that direct usage of this mass flux as the top boundary condition for the transient-diffusive
 *      vapor transport solver is not feasible, both theoretically and practically. This is primarily due to
 *      ql's turbulent nature as a latent heat flux, while the equation lacks a convection term.
 * -#   It calculates water vapor diffusivity and mass tranfer coefficient.
 * -#   When selecting the Explicit method, sub time steps are computed to ensure a stable solution. \n
 *      The method then integrates compDensityProfile to complete the full SNOWPACK time step, typically set at 15 minutes. \n
 *      Finally, mass is explicitly added or subtracted for each element, while adhering to specific constraints. \n
 * -#   For more information, please check the main reference as: DOI={10.3389/feart.2020.00249}
 * @author Mahdi Jafari
 * @param Xdata
 * @param ql Latent heat flux (W m-2)
 * @param Sdata
 * @param Mdata
 */
void VapourTransport::LayerToLayer(const CurrentMeteo& Mdata, SnowStation& Xdata, SurfaceFluxes& Sdata, double& ql)
{
	// First consider surface sublimation
	compSurfaceSublimation(Mdata, ql, Xdata, Sdata);

	const size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	std::vector<double> deltaM(nE, 0.);			// calculate the limited layer mass change

	if (!enable_vapour_transport) {
		// Only deal with the remaining ql (i.e., latent heat exchange at the surface)
		const double topFlux = -ql / Constants::lh_sublimation;										//top layer flux (kg m-2 s-1)
		const double dM = std::max(-EMS[nE-1].theta[ICE] * (Constants::density_ice * EMS[nE-1].L), -(topFlux * sn_dt));
		// Correct latent heat flux, which should become 0. at this point. HACK: note that if we cannot satisfy the ql at this point, we overestimated the latent heat from soil.
		// We will not get mass from deeper layers, as to do that, one should work with enable_vapour_transport == true.
		ql -= dM / sn_dt * Constants::lh_sublimation;
		deltaM[nE-1] += dM;
	} else {
		ql=0;

		size_t e = nE;
		std::vector<double> totalMassChange(nE, 0.);	// store the total mass change
		std::vector<double> oldVaporDenNode(nN, 0.);	// old water vapor density for node

		std::vector<double> factor_(nE, 1.);		// this is for source term in vapor transport equation
		for (size_t i = 0; i < Xdata.SoilNode; i++) {
			factor_[i] = 0.;
		}

		std::vector<double> D_(nE, 0.);
		for (size_t i = 0; i <= nE-1; i++) {
			double theta_air = std::max(EMS[i].theta[AIR], 0.0);
			double tortuosity = pow(theta_air, 7./3.) / pow(1-EMS[i].theta[SOIL], 2.);
			double D_vapSoil = tortuosity * theta_air * Constants::diffusion_coefficient_in_air;

			// based on Foslien (1994)
			double Dsnow = EMS[i].theta[ICE] * theta_air * Constants::diffusion_coefficient_in_air
						   + theta_air * Constants::diffusion_coefficient_in_air
						   / (EMS[i].theta[ICE] * Constants::conductivity_air
						   	  / Constants::conductivity_ice + EMS[i].theta[ICE]
							  * Constants::lh_sublimation * Constants::diffusion_coefficient_in_air
							  * dRhov_dT(EMS[i].Te) / Constants::conductivity_ice
							  + theta_air);
			D_[i] = factor_[i] * Dsnow + (1.0 - factor_[i]) * D_vapSoil;
		}

		std::vector<double> hm_(nN, 0.); // mass transfer coefficient m s-1

		for (size_t i = 0; i < nN; i++) {
			double saturationDensity;
			saturationDensity = Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));
			hm_[i] = Constants::density_ice / saturationDensity / 9.7e9; // hm_experimental, Pirmin 2012, M_mm=as_all*hm_experimental*(rhov_sat-rhov)
		}

		std::vector<double> as_(nN, 0.); // the specific surface area m-1
		for (size_t i=0; i<nN; i++) {
			if (i == 0) {
				if (!useSoilLayers) {
					double rwTors_u = pow(EMS[i].theta[WATER] / EMS[i].theta[ICE] + 1., 1. / 3.);
					double apparentTheta = EMS[i].theta[ICE] + EMS[i].theta[WATER];
					as_[i] = 6.0 * apparentTheta / (0.001 * rwTors_u * EMS[i].ogs);
				} else {
					double rwTors_u = pow((EMS[i].theta[WATER] + EMS[i].theta[ICE]) / EMS[i].theta[SOIL] + 1., 1. / 3.);
					double apparentTheta = EMS[i].theta[SOIL] + EMS[i].theta[ICE] + EMS[i].theta[WATER];
					as_[i] = 6.0 * apparentTheta / (0.002 * rwTors_u * EMS[i].rg);
				}
			} else if (i > 0 && i < Xdata.SoilNode) {
				double rwTors_u = pow((EMS[i].theta[WATER] + EMS[i].theta[ICE]) / EMS[i].theta[SOIL] + 1., 1. / 3.);
				double rwTors_d = pow((EMS[i-1].theta[WATER] + EMS[i-1].theta[ICE]) / EMS[i-1].theta[SOIL] + 1., 1. / 3.);
				double apparentTheta = 0.5 * (EMS[i].theta[SOIL] + EMS[i].theta[ICE] + EMS[i].theta[WATER]) + 0.5 * (EMS[i-1].theta[SOIL] + EMS[i-1].theta[ICE] + EMS[i-1].theta[WATER]);
				as_[i] = 6.0 * apparentTheta / (0.5 * 0.002 * rwTors_d * EMS[i-1].rg + 0.5 * 0.002 * rwTors_u * EMS[i].rg);
			} else if (i == Xdata.SoilNode && Xdata.SoilNode == nN-1) {
				double rwTors_d = pow((EMS[i-1].theta[WATER] + EMS[i-1].theta[ICE]) / EMS[i-1].theta[SOIL] + 1., 1./3.);
				double apparentTheta = EMS[i-1].theta[SOIL] + EMS[i-1].theta[ICE] + EMS[i-1].theta[WATER];
				as_[i] = 6.0 * apparentTheta / (0.002 * rwTors_d * EMS[i-1].rg);
			} else if (i == Xdata.SoilNode && Xdata.SoilNode < nN-1) {
				double rwTori_u = pow(EMS[i].theta[WATER] / EMS[i].theta[ICE] + 1., 1. / 3.);
				double rwTors_d = pow((EMS[i-1].theta[WATER] + EMS[i-1].theta[ICE]) / EMS[i-1].theta[SOIL] + 1., 1. / 3.);
				double apparentTheta = 0.5 * (EMS[i].theta[ICE] + EMS[i].theta[WATER]) + 0.5 * (EMS[i-1].theta[SOIL] + EMS[i-1].theta[ICE] + EMS[i-1].theta[WATER]);
				as_[i] = 6.0 * apparentTheta / (0.5 * 0.002 * rwTors_d * EMS[i-1].rg + 0.5 * 0.001 * rwTori_u * EMS[i].ogs);
			} else if (i > Xdata.SoilNode && i < nN-1) {
				double rwTori_u = pow(EMS[i].theta[WATER] / EMS[i].theta[ICE] + 1., 1. / 3.);
				double rwTori_d = pow(EMS[i-1].theta[WATER] / EMS[i-1].theta[ICE] + 1., 1. / 3.);
				double apparentTheta = 0.5 * (EMS[i].theta[ICE] + EMS[i].theta[WATER]) + 0.5 * (EMS[i-1].theta[ICE] + EMS[i-1].theta[WATER]);
				as_[i] = 6.0 * apparentTheta / (0.5 * 0.001 * rwTori_d * EMS[i-1].ogs + 0.5 * 0.001 * rwTori_u * EMS[i].ogs);
			} else { // i==nN-1
				double rwTori_d = pow(EMS[i-1].theta[WATER] / EMS[i-1].theta[ICE] + 1., 1. / 3.);
				double apparentTheta = EMS[i-1].theta[ICE] + EMS[i-1].theta[WATER];
				as_[i] = 6.0 * apparentTheta / (0.001 * rwTori_d * EMS[i-1].ogs);
			}
		}

		double min_dt = sn_dt; // first guess for the required minimum time step is the SNOWPACK time step
		for (size_t i=Xdata.SoilNode; i<nE; i++) {
			double saturationDensity = Atmosphere::waterVaporDensity(EMS[i].Te, Atmosphere::vaporSaturationPressure(EMS[i].Te));
			double diffVaporVelocity = std::abs(-D_[i] * (NDS[i+1].rhov-NDS[i].rhov) / EMS[i].L / saturationDensity);
			if (diffVaporVelocity!=0.) {
				double dt = EMS[i].L / diffVaporVelocity;
				min_dt = std::min(min_dt, dt);
			}
		}

		timeStep = (waterVaporTransport_timeStepAdjust) ? std::min(min_dt, waterVaporTransport_timeStep) : waterVaporTransport_timeStep;
		int nTime = int(sn_dt/timeStep) + 1;
		double time = 0.;
		for (size_t l = 0; static_cast<int>(l) <= nTime; l++) {
			time=time+timeStep;
			if (time >= sn_dt) {
				timeStep = sn_dt - (time - timeStep);
				time = sn_dt;
			}

			if (!compDensityProfile(Mdata, Xdata, hm_, as_, D_, oldVaporDenNode)) break;

			for (size_t i = 0; i <= nE-1; i++) {
				double saturationVaporUp = Atmosphere::waterVaporDensity(NDS[i+1].T, Atmosphere::vaporSaturationPressure(NDS[i+1].T));
				double saturationVaporDown = Atmosphere::waterVaporDensity(NDS[i].T, Atmosphere::vaporSaturationPressure(NDS[i].T));
				double diffRhov_hm_as_Up = (f * NDS[i+1].rhov + (1 - f) * oldVaporDenNode[i+1] - saturationVaporUp) * hm_[i+1] * as_[i+1];
				double diffRhov_hm_as_Down = (f * NDS[i].rhov + (1 - f) * oldVaporDenNode[i] - saturationVaporDown) * hm_[i] * as_[i];
				totalMassChange[i] = (0.5 * diffRhov_hm_as_Down + 0.5 * diffRhov_hm_as_Up) * timeStep * EMS[i].L; //total mass change, (kg m-2 )
			}

			e = nE;
			// consider the mass change due to vapour transport in snow/soil
			while (e-- > 0) {
				const double massPhaseChange = totalMassChange[e]+deltaM[e];

				double dM = 0.;	// mass change induced by vapor flux (kg m-2)

				// Now, the mass change is limited by:
				// - we cannot remove more WATER and ICE than available
				// - we cannot add more WATER and ICE than pore space available
				if ( EMS[e].theta[SOIL] < Constants::eps ) {	// there is no soil in element to keep element not to merge
					dM = std::max(-((EMS[e].theta[WATER] - EMS[e].VG.theta_r * (1. + Constants::eps)) * Constants::density_water * EMS[e].L + (EMS[e].theta[ICE] - Snowpack::min_ice_content) * Constants::density_ice * EMS[e].L),
								  std::min((EMS[e].theta[AIR] * Constants::density_ice * EMS[e].L), massPhaseChange)
						 		 ); // mass change due to difference in water vapor flux (kg m-2), at most can fill the pore space.
				} else {
					dM = std::max(-((EMS[e].theta[WATER] - EMS[e].VG.theta_r * (1. + Constants::eps)) * Constants::density_water * EMS[e].L + EMS[e].theta[ICE] * Constants::density_ice * EMS[e].L),
								  std::min((EMS[e].theta[AIR] * Constants::density_ice * EMS[e].L), massPhaseChange)
								 ); // mass change due to difference in water vapor flux (kg m-2), at most can fill the pore space.
				}


				// If there is no pore space, or, in fact, only so much pore space to accomodate the larger volume occupied by ice when all water freezes,
				// we inhibit vapour flux. This is necessary to maintain saturated conditions when present, and this is in turn necessary for the stability in the Richards equation solver.
				if (EMS[e].theta[AIR] < EMS[e].theta[WATER] * (Constants::density_water / Constants::density_ice - 1.) + Constants::eps) {
					dM = 0.;
				}

				deltaM[e] = dM;
			}

			if (time==sn_dt) break;
		}

		for (size_t i = 0; i < nE; i++) {
			EMS[i].vapTrans_fluxDiff = -D_[i] * (NDS[i+1].rhov-NDS[i].rhov) / EMS[i].L;
		}

		double dHoar = 0.;
		for (e = 0; e < nE; e++) {
			EMS[e].Qmm = 0.0;

			if (deltaM[e] < 0.) {
				// Mass loss: apply mass change first to water, then to ice, based on energy considerations
				// We can only do this partitioning here in this "simple" way, without checking if the mass is available, because we already limited dM above, based on available ICE + WATER.
				const double dTh_water = std::max((EMS[e].VG.theta_r * (1. + Constants::eps) - EMS[e].theta[WATER]),
												  deltaM[e] / (Constants::density_water * EMS[e].L));
				const double dTh_ice = ( deltaM[e] - (dTh_water * Constants::density_water * EMS[e].L) ) / (Constants::density_ice * EMS[e].L);
				EMS[e].theta[WATER] += dTh_water;
				EMS[e].theta[ICE] += dTh_ice;

				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dTh_water * Constants::density_water * EMS[e].L;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dTh_ice * Constants::density_ice * EMS[e].L;
				EMS[e].M += dTh_water * Constants::density_water * EMS[e].L+dTh_ice * Constants::density_ice * EMS[e].L;
				assert(EMS[e].M >= (-Constants::eps2)); // mass must be positive

				EMS[e].Qmm += (dTh_water * Constants::density_water * Constants::lh_vaporization
							   + dTh_ice * Constants::density_ice * Constants::lh_sublimation
							  ) / sn_dt; // [w/m^3]

				// If present at surface, surface hoar is sublimated away
				if (e == nE-1 && deltaM[e]<0) {
					dHoar = std::max(-NDS[nN-1].hoar, deltaM[e]);
				}
			} else {  // Mass gain: add water in case temperature at or above melting point, ice otherwise
				if (EMS[e].Te >= EMS[e].meltfreeze_tk) {
					EMS[e].theta[WATER] += deltaM[e] / (Constants::density_water * EMS[e].L);
					EMS[e].Qmm += (deltaM[e]*Constants::lh_vaporization)/sn_dt/EMS[e].L;	// [w/m^3]
					Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += deltaM[e];
				} else {
					EMS[e].theta[ICE] += deltaM[e] / (Constants::density_ice * EMS[e].L);
					EMS[e].Qmm += (deltaM[e]*Constants::lh_sublimation)/sn_dt/EMS[e].L;	// [w/m^3]
					Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += deltaM[e];
				}
				EMS[e].M += deltaM[e];
				assert(EMS[e].M >= (-Constants::eps2)); // mass must be positive
			}


			EMS[e].theta[AIR] = std::max(1. - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL], 0.);
			if (std::fabs(EMS[e].theta[AIR]) < 1.e-15) {
				EMS[e].theta[AIR] = 0;
			}
			EMS[e].updDensity();
			assert(EMS[e].Rho > 0 || EMS[e].Rho == IOUtils::nodata); // density must be positive
			if (!(EMS[e].Rho > Constants::eps && EMS[e].theta[AIR] >= 0. && EMS[e].theta[WATER] <= 1. + Constants::eps && EMS[e].theta[ICE] <= 1. + Constants::eps)) {
					prn_msg(__FILE__, __LINE__, "err", Date(),
						"Volume contents: e=%d nE=%d rho=%lf ice=%lf wat=%lf wat_pref=%lf soil=%lf air=%le", e, nE, EMS[e].Rho, EMS[e].theta[ICE],
							EMS[e].theta[WATER], EMS[e].theta[WATER_PREF], EMS[e].theta[SOIL], EMS[e].theta[AIR]);
					throw IOException("Cannot evaluate mass balance in vapour transport LayerToLayer routine", AT);
			}

			// some useful output in case of vapor transport
			double sVaporDown = Atmosphere::waterVaporDensity(NDS[e].T, Atmosphere::vaporSaturationPressure(NDS[e].T));
			double sVaporUp = Atmosphere::waterVaporDensity(NDS[e+1].T, Atmosphere::vaporSaturationPressure(NDS[e+1].T));
			EMS[e].vapTrans_underSaturationDegree = (0.5*(NDS[e].rhov-sVaporDown)+0.5*(NDS[e+1].rhov-sVaporUp))/(0.5*sVaporDown+0.5*sVaporUp);
			EMS[e].vapTrans_cumulativeDenChange += deltaM[e]/EMS[e].L;
			EMS[e].vapTrans_snowDenChangeRate = deltaM[e]/EMS[e].L/sn_dt;
		}

		Sdata.hoar += dHoar;
		NDS[nN-1].hoar += dHoar;
		if (NDS[nN-1].hoar < 0.) {
			NDS[nN-1].hoar = 0.;
		}

	}
}

/**
 * @brief Calculate the surface sublimation / deposition (i.e., only gas-solid). \n
 * The fraction of the latent heat flux ql that has not been used so far will be used for
 * sublimation/deposition. If positive (and above a certain cutoff level) then there
 * is a possibility that surface hoar crystal have grown. Of course, if negative
 * then we are also loosing mass from the surface.\n
 * This function additionally takes care of surface hoar formation and destruction.
 * Note that surface hoar is a nodal property, altough the corresponding mass is carried
 * by the underlying element.
 * @param *Mdata
 * @param ql Latent heat flux (W m-2)
 * @param *Xdata
 * @param *Sdata
 */
void VapourTransport::compSurfaceSublimation(const CurrentMeteo& Mdata, double& ql, SnowStation& Xdata, SurfaceFluxes& Sdata)
{
	double dM, M;		// Length and mass changes, and Initial mass and volumetric content (water or ice)
	double dHoar = 0.;	// Actual change in hoar mass
	double cH_old;		// Temporary variable to hold height of snow

	const size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	const double Tss = NDS[nE].T; // Surface Temperature

	/*
	 * If ql > 0:
	 * Surface hoar is formed when surface temperature is below freezing.
	 * If no surface hoar can be formed, ql is kept and is used as boundary condition
	 * when calculating vapour flux.
	 * If there are elements and ql < 0:
	 * If ql is large enough to remove full surface elements, remove them.
	 * left over ql is used as boundary condition when calculating vapour flux.
	 *
	 * In both cases: add/subtract mass to MS_SUBLIMATION
	 */
	if (ql > Constants::eps2) { // Add Mass
		const double meltfreeze_tk = (Xdata.getNumberOfElements()>0)? Xdata.Edata[Xdata.getNumberOfElements()-1].meltfreeze_tk : Constants::meltfreeze_tk;
		if (Tss < meltfreeze_tk) { // Add Ice
			dM = ql*sn_dt/Constants::lh_sublimation;
			// If rh is very close to 1, vw too high or ta too high, surface hoar is destroyed and should not be formed
			if (!((Mdata.rh > hoar_thresh_rh) || (Mdata.vw > hoar_thresh_vw) || (Mdata.ta >= IOUtils::C_TO_K(hoar_thresh_ta)))) {
				// Under these conditions, form surface hoar
				ql = 0.;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
				dHoar = dM;

				// In this case adjust properties of element, keeping snow density constant
				const double L_top = EMS[nE-1].L;
				const double theta_i0 = EMS[nE-1].theta[ICE];
				double dL = dM/(EMS[nE-1].Rho); // length change
				if (nE == Xdata.SoilNode) {
					dL = 0.;
					dM = std::min(dM,EMS[nE-1].theta[AIR]*(Constants::density_ice*EMS[nE-1].L));
				}
				NDS[nE].z += dL + NDS[nE].u; NDS[nE].u = 0.0;
				EMS[nE-1].L0 = EMS[nE-1].L = L_top + dL;
				EMS[nE-1].E = EMS[nE-1].Eps = EMS[nE-1].dEps = EMS[nE-1].Eps_e = EMS[nE-1].Eps_v = EMS[nE-1].S = 0.0;
				EMS[nE-1].theta[ICE] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[ICE] += dM/(Constants::density_ice*EMS[nE-1].L);
				EMS[nE-1].theta[ICE] = std::max(0., std::min(1., EMS[nE-1].theta[ICE]));
				EMS[nE-1].theta[WATER] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[WATER] = std::max(0., std::min(1., EMS[nE-1].theta[WATER]));
				EMS[nE-1].theta[WATER_PREF] *= L_top/EMS[nE-1].L;
				EMS[nE-1].theta[WATER_PREF] = std::max(0., std::min(1., EMS[nE-1].theta[WATER_PREF]));

				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[nE-1].conc[ICE][ii] *= L_top*theta_i0/(EMS[nE-1].theta[ICE]*EMS[nE-1].L);
				}

				EMS[nE-1].M += dM;
				assert(EMS[nE-1].M >= (-Constants::eps2)); //mass must be positive

				// Update remaining volumetric contents and density
				EMS[nE-1].theta[AIR] = std::max(0., 1.0 - EMS[nE-1].theta[WATER] - EMS[nE-1].theta[WATER_PREF] - EMS[nE-1].theta[ICE] - EMS[nE-1].theta[SOIL]);
				EMS[nE-1].updDensity();
			}
		}
	} else if ((ql < (-Constants::eps2)) && (nE > 0)) {
		// If ql < 0, SUBLIMATE mass off
		std::vector<double> M_Solutes(Xdata.number_of_solutes, 0.); // Mass of solutes from disappearing phases
		size_t e = nE;
		while ((e > 0) && (ql < (-Constants::eps2))) {	// While energy is available
			e--;
			/*
			* Determine the amount of potential sublimation and collect some variables
			* that will be continuously used: L0 and M
			*/
			const double L0 = EMS[e].L;
			const double theta_i0 = EMS[e].theta[ICE];

			M = theta_i0 * Constants::density_ice * L0;
			dM = ql * sn_dt / Constants::lh_sublimation;

			if (-dM > M) {
				// Only if mass change is sufficient to remove the full element
				dM = -M;
				// Add solutes to Storage
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					M_Solutes[ii] += EMS[e].conc[ICE][ii]*theta_i0*L0;
				}
				EMS[e].theta[ICE] = 0.;

				EMS[e].M += dM;
				Sdata.mass[SurfaceFluxes::MS_SUBLIMATION] += dM;
				ql -= dM*Constants::lh_sublimation/sn_dt;	// Update the energy used

				// If present at surface, surface hoar is sublimated away
				if (e == nE-1) {
					dHoar = std::max(-NDS[nN-1].hoar, dM);
				}

				// Update remaining volumetric contents and density
				EMS[e].theta[AIR] = std::max(0., 1.0 - EMS[e].theta[WATER] - EMS[e].theta[WATER_PREF] - EMS[e].theta[ICE] - EMS[e].theta[SOIL]);
				EMS[e].updDensity();
				// Merge the element if it is a snow layer. This will take care of possible left over liquid water (will be put one layer down)
				// Keep layer if it is a soil layer inside the snowpack (for example with snow farming)
				if (e >= Xdata.SoilNode) {
					if (EMS[e].theta[SOIL] < Constants::eps) {
						if (e > 0) SnowStation::mergeElements(EMS[e-1], EMS[e], false, true);
						// Now reduce the number of elements by one.
						nE--;
					}
					//In case e==Xdata.SoilNode, we removed the last snow element and we should break out of the loop.
					if (e == Xdata.SoilNode) break;
				}
			} else {
				// Not enough energy anymore to remove complete element, so we should break out of the loop.
				break;
			}

			//check that thetas and densities are consistent
			assert(EMS[e].theta[SOIL] >= (-Constants::eps2) && EMS[e].theta[SOIL] <= (1.+Constants::eps2));
			assert(EMS[e].theta[ICE] >= (-Constants::eps2) && EMS[e].theta[ICE]<=(1.+Constants::eps2));
			assert(EMS[e].theta[WATER] >= (-Constants::eps2) && EMS[e].theta[WATER]<=(1.+Constants::eps2));
			assert(EMS[e].theta[WATER_PREF] >= (-Constants::eps2) && EMS[e].theta[WATER_PREF]<=(1.+Constants::eps2));
			assert(EMS[e].theta[AIR] >= (-Constants::eps2) && EMS[e].theta[AIR]<=(1.+Constants::eps2));
			assert(EMS[e].Rho >= (-Constants::eps2) || EMS[e].Rho==IOUtils::nodata); //we want positive density
		}

		// Now take care of left over solute mass.
		if (nE == Xdata.SoilNode) { // Add Solute Mass to Runoff TODO HACK CHECK
			for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
				Sdata.load[ii] += M_Solutes[ii]/S_TO_H(sn_dt);
			}
		} else { // Add Solute Mass to Element below
			if (EMS[e].theta[WATER] > 0.) {
				for(size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[WATER][ii] += M_Solutes[ii]/EMS[e].theta[WATER]/EMS[e].L;
				}
			} else if (EMS[e].theta[ICE] > 0.) {
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[ICE][ii] += M_Solutes[ii]/EMS[e].theta[ICE]/EMS[e].L;
				}
			} else {
				for (size_t ii = 0; ii < Xdata.number_of_solutes; ii++) {
					EMS[e].conc[SOIL][ii] += M_Solutes[ii]/EMS[e].theta[SOIL]/EMS[e].L;
				}
			}
		}
		Xdata.reduceNumberOfElements(nE);
	}

	// HACK: this code is under verification. The comment reads "surface hoar *is* destroyed, but the next line says surface hoar *may be* destroyed, depending on the sign of the latent heat flux.
	// If the code is correct, we can delete this part, if the comment is correct, we should modify the code to read: hoar = -NDS[nE].hoar;
	// Check for surface hoar destruction or formation (once upon a time ml_sn_SurfaceHoar)
	/*if ((Mdata.rh > hoar_thresh_rh) || (Mdata.vw > hoar_thresh_vw) || (Mdata.ta >= IOUtils::C_TO_K(hoar_thresh_ta))) {
		//if rh is very close to 1, vw too high or ta too high, surface hoar is destroyed
		hoar = std::min(hoar, 0.);
	}*/

	Sdata.hoar += dHoar;
	NDS[nN-1].hoar += dHoar;
	if (NDS[nN-1].hoar < 0.) {
		NDS[nN-1].hoar = 0.;
	}

	// Surface hoar cannot exist when the top element is wet
	if (nE > 0) {
		const double theta_r = ((iwatertransportmodel_snow == RICHARDSEQUATION && nE-1>=Xdata.SoilNode) || (iwatertransportmodel_soil==RICHARDSEQUATION && nE-1<Xdata.SoilNode)) ? (PhaseChange::RE_theta_r) : (PhaseChange::theta_r);
		if (Xdata.Edata[nE-1].theta[WATER] > theta_r) {
			NDS[nE].hoar = 0.;
		}
	}

	// At the end also update the overall height
	cH_old = Xdata.cH;
	Xdata.cH = NDS[Xdata.getNumberOfNodes()-1].z + NDS[Xdata.getNumberOfNodes()-1].u;
	if (Xdata.mH!=Constants::undefined) Xdata.mH -= std::min(Xdata.mH - Xdata.Ground, (cH_old - Xdata.cH));	// TODO/HACK: why is this correction for Xdata.mH necessary?
}

/**
 * @brief This function is the solver for discretized transient-diffusive vapor tranport equation.
 * NOTES:
 * -#   Note, for the case of only snow (no soil), bottomDirichletBCtype is set to Drichlet ans Neumann does not make sense \n
 * -#   The system of equations forms a tridiagonal sparse matrix for which the sparse solvers from the Eigen C++ library are used. \n
 *      Here, we used quite well stabel solver as BiCGSTAB. Feel free to use other solvers by looking at Eigen documentaion.
 * -#   When selecting the Explicit method, sub time steps are computed to ensure a stable solution. \n
 *      The method then integrates compDensityProfile to complete the full SNOWPACK time step, typically set at 15 minutes. \n
 *      Finally, mass is explicitly added or subtracted for each element, while adhering to specific constraints. \n
 * @author Mahdi Jafari
 * @param Xdata
 * @param Mdata
 * @param D_el, water vapor diffusivity (m2 s-1)
 * @param hm_, mass tranfer coefficient (m s-1)
 * @param as_, specific surface area (m-1)
 * @param oldVaporDenNode, old vapor denisty stored for nodes
 * @return a flag to check if the solution has converged
 */
bool VapourTransport::compDensityProfile(const CurrentMeteo& Mdata, SnowStation& Xdata,
										 std::vector<double>& hm_,
										 std::vector<double>& as_,
										 const std::vector<double>& D_el,
										 std::vector<double>& oldVaporDenNode)
{
	const bool bottomDirichletBCtype = (Xdata.SoilNode == 0 && variant != "SEAICE") ? (true) : (false);

	const size_t nN = Xdata.getNumberOfNodes();
	size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;
	const size_t nX = nN; // number of unknowns

	BiCGSTAB<SparseMatrix<double> > solver; // Built-in iterative solver

	SparseMatrix<double,RowMajor> A(nX,nX);
	std::vector<Trip> tripletList(nX);
	VectorXd b(nX);
	VectorXd xx(nX);

	// grid
	std::vector<double> z(nN,0.);
	for (size_t i = 0; i < nN; i++) {
		z[i] = NDS[i].z;
	}

	// initial values
	std::vector<double> D_(nN, Constants::diffusion_coefficient_in_air);
	for(size_t i=0; i<nN; i++) {
		if (i == 0) {
			D_[i]=D_el[i];
		} else if (i==nN-1) {
			D_[i]=0.5*D_el[i-1]+0.5*Constants::diffusion_coefficient_in_air;
		} else {
			D_[i]=D_el[i];
		}
	}

	// initial values
	// if the effective water vapor diffusivity in air is defined eps_=thata_a,
	// otherwise for the effective water vapor diffusivity in snow eps_=1.
	// for now Deff,s is available so eps_=1
	std::vector<double> eps_(nN, 1.0);
	for(size_t i=0; i<nN; i++)
	{
		if(i==0)
		{
			eps_[i]=std::max(EMS[i].theta[AIR],1.0);
		}
		else if(i==nN-1)
		{
			eps_[i]=0.5*std::max(EMS[i-1].theta[AIR],1.0)+0.5;
		}
		else
		{
			eps_[i]=std::max(EMS[i].theta[AIR],1.0);
		}
	}

	double error_max = 0;
	do
	{
		error_max = 0;

		// The lower B.C.
		if(bottomDirichletBCtype){
			double elementSaturationVaporDensity=Atmosphere::waterVaporDensity(NDS[0].T, Atmosphere::vaporSaturationPressure(NDS[0].T));
			NDS[0].rhov=elementSaturationVaporDensity;
		}

		// Diffusion equation
		double v_ij;
		A.setZero();
		double dz_u, dz_d;
		double saturationDensity;
		double eps_n;

		for(size_t k=0; k<=nN-1; k++) {
			saturationDensity = Atmosphere::waterVaporDensity(NDS[k].T, Atmosphere::vaporSaturationPressure(NDS[k].T));

			if (k != 0 && k != nN-1) {
				dz_u = z[k+1]-z[k];
				dz_d = z[k]-z[k-1];
				eps_n = (0.5 * eps_[k] + 0.5 * eps_[k-1]);

				b[k] = eps_n*NDS[k].rhov/timeStep+hm_[k]*as_[k]*(saturationDensity-(1-f)*NDS[k].rhov);
				b[k] += -NDS[k].rhov*((1-f)*2.0*eps_[k]*D_[k]/dz_u/(dz_u+dz_d)+(1-f)*2.0*eps_[k-1]*D_[k-1]/dz_d/(dz_u+dz_d));
				b[k] += -NDS[k+1].rhov*(1-f)*-2.0*eps_[k]*D_[k]/dz_u/(dz_u+dz_d);
				b[k] += -NDS[k-1].rhov*(1-f)*-2.0*eps_[k-1]*D_[k-1]/dz_d/(dz_u+dz_d);

				v_ij = f * 2.0 * eps_[k] * D_[k] / dz_u / (dz_u + dz_d) + f * 2.0 * eps_[k-1] * D_[k-1] / dz_d / (dz_u + dz_d) + eps_n * 1.0 / timeStep + hm_[k] * as_[k] * f;
				tripletList.push_back(Trip(static_cast<int>(k), static_cast<int>(k), v_ij));		// Set up the matrix diagonal

				v_ij = f * -2.0 * eps_[k] * D_[k] / dz_u / (dz_u + dz_d);
				tripletList.push_back(Trip(static_cast<int>(k), static_cast<int>(k) + 1, v_ij));	// Set up the matrix upper diagonals, k+1

				v_ij = f * -2.0 * eps_[k-1] * D_[k-1] / dz_d / (dz_u + dz_d);
				tripletList.push_back(Trip(static_cast<int>(k), static_cast<int>(k) - 1, v_ij));	// Set up the matrix lower diagonals, k-1
			} if (k == nN-1) {
				// Normal top B.C. assuming satuarion condition for the uppermost node of snowpack
				b[k] = saturationDensity;
				v_ij = 1.0;
				tripletList.push_back(Trip(static_cast<int>(k), static_cast<int>(k), v_ij));		// Set up the matrix diagonal
			} if (k == 0) {
				if (bottomDirichletBCtype) {
					b[k] = saturationDensity;  // NDS[k].rhov;
					v_ij = 1.0;
					tripletList.push_back(Trip(static_cast<int>(k), static_cast<int>(k), v_ij));	// Set up the matrix diagonal
				} else {
					b[k] = 0.0;
					v_ij = -1.0;
					tripletList.push_back(Trip(static_cast<int>(k), static_cast<int>(k), v_ij));	// Set up the matrix diagonal

					v_ij = 1.0;
					tripletList.push_back(Trip(static_cast<int>(k), static_cast<int>(k) + 1, v_ij));// Set up the matrix upper diagonals, k+1
				}
			}
		}

		A.setFromTriplets(tripletList.begin(), tripletList.end());
		tripletList.clear();
		A.makeCompressed();

		solver.compute(A);
		if (solver.info() != Success) {
			std::ostringstream err_msg;
			err_msg << "Error computing 'A' with Eigen: " << Mdata.date << solver.info();
			throw mio::IOException(err_msg.str(), AT);
		}

		// Solve the equation
		xx = solver.solve(b);
		if (solver.info() != Success) {
			std::ostringstream err_msg;
			err_msg << "Error solving 'b' with Eigen: " << Mdata.date << solver.info();
			throw mio::IOException(err_msg.str(), AT);
		}

		for (size_t k = 0; k <= nN-1; k++) {
			oldVaporDenNode[k]=NDS[k].rhov;
			NDS[k].rhov=xx(k);
			double error = std::abs(NDS[k].rhov-oldVaporDenNode[k]);
			if(NDS[k].rhov<0) {
				std::ostringstream err_msg;
				err_msg << "Error, rhov is below zero (" << NDS[k].rhov << "). Can not proceed.";
				throw mio::IOException(err_msg.str(), AT);
			}
			error_max = std::max(error_max, error);
			saturationDensity = Atmosphere::waterVaporDensity(NDS[k].T, Atmosphere::vaporSaturationPressure(NDS[k].T));
		}

		break;
	} while (error_max > 1.e-6);

	for (size_t e = 0; e < nE; e++) {
		EMS[e].rhov = (NDS[e].rhov + NDS[e+1].rhov) / 2.0;
	}

	return true;
}

/**
 * @brief Calculate the derivative of vapor saturation pressure of a flat ice surface with respect to temperature.\n
 * Refer to MeteoIO::Atmosphere::vaporSaturationPressure() for the non-derived method
 * in the MeteoIO library.
 * @param Tem Temperature in Kelvin
 * @return Derivative of vapor saturation pressure with respect to temperature
 */
double VapourTransport::dRhov_dT(const double Tem)
{
	/* Use the constants of a flat ice surface.
	 * See Murray, F. W., "On the computation of saturation vapor pressure", 1966, J. Appl. Meteor., 6, 203–204,
	 * doi: 10.1175/1520-0450(1967)006<0203:OTCOSV>2.0.CO;2.
	 */
	const double sat_vapor_pressure_ice_a = 21.8745584;
	const double sat_vapor_pressure_ice_b = 7.66;

	const double dRhov_dT = (mio::Cst::water_molecular_mass * mio::Cst::p_water_triple_pt / mio::Cst::gaz_constant / Tem)
		* exp(sat_vapor_pressure_ice_a * (Tem - mio::Cst::t_water_triple_pt) / (Tem - sat_vapor_pressure_ice_b))
		* (-1. / Tem + (sat_vapor_pressure_ice_a * (Tem - sat_vapor_pressure_ice_b) - sat_vapor_pressure_ice_a * (Tem - mio::Cst::t_water_triple_pt)) / (Tem - sat_vapor_pressure_ice_b) / (Tem - sat_vapor_pressure_ice_b));

	return dRhov_dT;
}
/*
 * End of VapourTransport.cc
 */
