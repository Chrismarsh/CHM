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

#include <snowpack/Stability.h>

#include <assert.h>

using namespace mio;
using namespace std;

/************************************************************
 * static section                                           *
 ************************************************************/

const double Stability::psi_ref = 38.0; ///< Reference slope angle
const double Stability::max_stability = 6.0; ///< Upper stability limit

///Minimum slab thickness for natural and deformation stability index (m)
const double Stability::minimum_slab = 0.1;

///The first GROUND_ROUGH m of snow will not be unstable due to ground roughness
const double Stability::ground_rough = 0.2;

///MIN_DEPTH_SSI m of snow must be left after discarding penetration depth
const double Stability::min_depth_ssi = 0.1;

///Skiers will not trigger failures SKIER_DEPTH m below penetration depth
const double Stability::skier_depth = 1.0;

///Minimum thickness for a supporting melt-freeze crust (perp to slope, in m)
const double Stability::min_thick_crust = 0.03;

///Maximum number of structural instabilities looked at ("lemons")
const size_t Stability::nmax_lemon = 2;

///Defines regression model for surface hoar shear strength
const int Stability::sh_mod = 2;

/**
 * @brief Defines classification scheme for snow profiles
 * - 0: Based on Skier Stability Index (SSI) thresholds (3 classes)
 * - 1: Based on re-analysis by Schweizer/Bellaire of SSI and SK38 (2006)
 * - 2: Nov 2007: re-analysis after recalibration of settling (see rev 250/251)
 * - 3: According to Schweizer and Wiesinger (5 classes)
 */
const int Stability::prof_classi = 2;

map<string, StabMemFn> Stability::mapHandHardness;
map<string, StabFnShearStrength> Stability::mapShearStrength;
const bool Stability::__init = Stability::initStaticData();

bool Stability::initStaticData()
{
	mapHandHardness["MONTI"]    = &Stability::setHandHardnessMONTI;
	mapHandHardness["BELLAIRE"]  = &Stability::setHandHardnessBELLAIRE;
	mapHandHardness["ASARC"]    = &Stability::setHandHardnessASARC;

	mapShearStrength["DEFAULT"] = &Stability::setShearStrengthDEFAULT;
	mapShearStrength["NIED"]    = &Stability::setShearStrengthSTRENGTH_NIED;

	return true;
}

/************************************************************
 * non-static section                                       *
 ************************************************************/

Stability::Stability(const SnowpackConfig& cfg, const bool& i_classify_profile)
           : strength_model(), hardness_parameterization(), hoar_density_buried(IOUtils::nodata), plastic(false),
             classify_profile(i_classify_profile)
{
	cfg.getValue("STRENGTH_MODEL", "SnowpackAdvanced", strength_model);
	cfg.getValue("HARDNESS_PARAMETERIZATION", "SnowpackAdvanced", hardness_parameterization);

	const map<string, StabMemFn>::const_iterator it1 = mapHandHardness.find(hardness_parameterization);
	if (it1 == mapHandHardness.end()) throw InvalidArgumentException("Unknown hardness parameterization: "+hardness_parameterization, AT);

	const map<string, StabFnShearStrength>::const_iterator it2 = mapShearStrength.find(strength_model);
	if (it2 == mapShearStrength.end()) throw InvalidArgumentException("Unknown strength model: "+strength_model, AT);

	//To build a sandwich with a non-snow layer (plastic or wood chips) on top;
	//originally introduced for snow farming applications
	cfg.getValue("PLASTIC", "SnowpackAdvanced", plastic);

	// Density of BURIED surface hoar (kg m-3), default: 125./ Antarctica: 200.
	cfg.getValue("HOAR_DENSITY_BURIED", "SnowpackAdvanced", hoar_density_buried);
}


/**
 * @brief Assign hardness to snow types according to density, Swiss version by Sascha Bellaire
 * @author Implemented by C. Fierz: Regression by Sascha Bellaire January 2005 (all types except MFcr).
 * The original Swiss regression has been modified for PP, RG and FC to get a better agreement
 * to observed hardness. If there is a new settlement formulation in future, and therefore a
 * better agreement to observed density,it must be checked whether the new or the original
 * regression is the better choice. (18 September 2005; Fierz / S. Bellaire)
 * Original regression values are added as comments where needed.
 * @param Edata
 * @return hand hardness index (1)
 */
double Stability::setHandHardnessBELLAIRE(const ElementData& Edata)
{
	const double gsz = 2.*Edata.rg;
	double hardness;

	if ( (Edata.mk%100) < 20 ) { // all types except MFcr (hardness 5)
		int F1, F2, F3; // grain shape
		typeToCode(&F1, &F2, &F3, Edata.type); // Decompose type in its constituents
		double A, B;
		switch(F1) {
			case 0: { // Graupel PPgp; introduced by Yamaguchi & Fierz, Feb 2004
				A = 0.0078;
				B = 0.0105;
				break;
			}
			case 1: { // Precipitation Particles PP; ori: A = 0.7927; B = 0.0038;
				A = 0.7927;
				B = 0.0036;
				break;
			}
			case 2: { // Decomposing and Fragmented precipitation particles DF
				A = 0.4967;
				B = 0.0074;
				break;
			}
			case 3: { // Rounded Grains RG; ori: A = 0.2027; B = 0.0092;
				A = 0.2027;
				B = 0.0072;
				break;
			}
			case 4: { // Faceted Crystals FC; ori: A = 0.3867; B = 0.0071;
				A = 0.3867;
				B = 0.0083;
				break;
			}
			case 5: { // Depth Hoar DH
				A = -0.00249;
				B = 0.0072;
				break;
			}
			case 6: { // Surface hoar SH; empirical: index 1 to 2 from hoar_density_buried to 250 kg m-3
				A = 1. - hoar_density_buried/(250. - hoar_density_buried);
				B = 1./(250. - hoar_density_buried);
				break;
			}
			case 7: { // Melt Forms MF
				A = 0.5852;
				B = 0.0056;
				break;
			}
			case 8: { // Ice layer IFil
				A = 6.;
				B = 0.;
				break;
			}
			case 9: { // Rounding faceted particles FCxr
				A = -0.5226;
				B = 0.0104;
				break;
			}
			default: {
				A = Constants::undefined;
				B = 0.;
				break;
			}
		}
		hardness = A + B*Edata.Rho;
		// Large surface hoar stays longer unstable! 1 dec 2007 (sb)
		if ((F1 == 6) && (gsz >= 5.)) {
			hardness = 1;
		} else if ((F1 == 6 ) && (gsz < 5.)) {
			hardness = MIN(hardness, 2.);
		}
	} else if (Edata.theta[ICE] <= 0.7) { // Melt-freeze crust MFcr
		if (Edata.theta[WATER] < 0.3 * Edata.res_wat_cont) {
			hardness = 5.;
		} else if (Edata.theta[WATER] < 0.6 * Edata.res_wat_cont) {
			hardness = 4.5;
		} else if (Edata.theta[WATER] < 0.85 * Edata.res_wat_cont) {
			hardness = 4.;
		} else {
			hardness = 3.;
		}
	} else { // Ice Formations IF
		hardness = 6.;
	}
	// Limit to range {1, 6}
	hardness = MAX(1., MIN(6., hardness));
	return hardness;
}

/**
 * @brief Compute hand hardness for a given grain type and density
 * All the information about hardness parameterizations for PP DF RG FC DH MF FCxf are published in
 * Monti et al. (in progress)
 * @author Fabiano Monti
 * @date 2012-06-27
 * @param F grain type
 * @param rho snow density
 * @return hand hardness index (1)
 */
double Stability::getHandHardnessMONTI(const int& F, const double& rho, const double& water_content)
{
	switch(F) {
		case 0: { // Graupel PPgp; introduced by Yamaguchi & Fierz, Feb 2004
			const double A = 0.0078;
			const double B = 0.0105;
			return (A + B*rho);
		}
		case 1: { // Precipitation Particles PP; obtained from median value for hand_hardness_1 (110 kg/m3) + standard dev (33.9397 kg/m3)
		          // if not the median value for hand_hardness_2 is 129.5 kg/m3 but it comes from only 6 observations;
			if ((rho >= 0.) && (rho < 143.9397))
				return 1.;
			else
				return 2.;
		}
		case 2: { // Decomposing and Fragmented precipitation particles DF
			if ((rho >= 0.) && (rho <= 214.2380))
				return 1.;
			else if ((rho > 214.2380) && (rho <= 268.2981))
				return 2.;
			else if ((rho > 268.2981) && (rho <= 387.4305))
				return 3.;
			else
				return 4.;
		}
		case 3: { // Rounded Grains RG
			if ((rho >= 0.) && (rho <= 189.2103))
				return 1.;
			else if ((rho > 189.2103) && (rho <= 277.8087))
				return 2.;
			else if ((rho > 277.8087) && (rho <= 368.4093))
				return 3.;
			else if ((rho > 368.4093) && (rho <= 442.4917))
				return 4.;
			else
				return 5.;
		}
		case 4: { // Faceted Crystals FC
			if ((rho >= 0.) && (rho <= 247.2748))
				return 1.;
			else if ((rho > 247.2748) && (rho <= 319.3549))
				return 2.;
			else if ((rho > 319.3549) && (rho <= 400.4450))
				return 3.;
			else if ((rho > 400.4450) && (rho <= 517.5751))
				return 4.;
			else
				return 5.;
		}
		case 5: { // Depth Hoar DH
			if ((rho >= 0.) && (rho <= 287.8198))
				return 1.;
			else if ((rho > 287.8198) && (rho <= 344.3826))
				return 2.;
			else
				return 3.;
		}
		case 6: { // Surface hoar SH; empirical: index 1 to 2 from hoar_density_buried to 250 kg m-3
			const double A = 1. - hoar_density_buried/(250. - hoar_density_buried);
			const double B = 1./(250. - hoar_density_buried);
			return (A + B*rho);
		}
		case 7: { // Melt Forms MF
			if (water_content < SnowStation::thresh_moist_snow) { //dry melt forms
				if ((rho >= 0.) && (rho <= 213.7375))
					return 1.;
				else if ((rho > 213.7375) && (rho <= 317.3527))
					return 2.;
				else if ((rho > 317.3527) && (rho <= 406.9522))
					return 3.;
				else if ((rho > 406.9522) && (rho <= 739.8220))
					return 4.;
				else
					return 5.;
			} else { //moist melt forms
				if ((rho >= 0.) && (rho <= 338.3760))
					return 1.;
				else if ((rho > 338.3760) && (rho <= 417.4638))
					return 2.;
				else if ((rho > 417.4638) && (rho <= 541.6018))
					return 3.;
				else if ((rho > 541.6018) && (rho <= 614.6830))
					return 4.;
				else
					return 5.;
			}
		}
		case 8: { // Ice layer IFil
			const double A = 6.;
			const double B = 0.;
			return (A + B*rho);
		}
		case 9: { // Rounding faceted particles FCxr
			if ((rho >= 0.) && (rho <= 259.7887))
				return 1.;
			else if ((rho > 259.7887) && (rho <= 326.8632))
				return 2.;
			else if ((rho > 326.8632) && (rho <= 396.9411))
				return 3.;
			else if ((rho > 396.9411) && (rho <= 484.5384))
				return 4.;
			else
				return 5.;
		}
		default: {
			std::stringstream ss;
			ss << "Error: grain type " << F << " is unknown!";
			throw IOException(ss.str(), AT);
		}
	}
}

/**
 * @brief Assign hardness to snow types according to density
 * Implementation according to Fabiano Monti's work, June 2012 (all types except MFcr).
 * @author Fabiano Monti
 * @date 2012-06-27
 * @param Edata
 * @return hand hardness index (1)
 */
double Stability::setHandHardnessMONTI(const ElementData& Edata)
{
	double hardness;

	if ( (Edata.mk%100) < 20 ) { // all types except MFcr (hardness 5)
		int F1, F2, F3; // grain shape
		typeToCode(&F1, &F2, &F3, Edata.type); // Decompose type in its constituents
		const double hardness_F1 = getHandHardnessMONTI(F1, Edata.Rho, Edata.theta[WATER]);
		const double hardness_F2 = getHandHardnessMONTI(F2, Edata.Rho, Edata.theta[WATER]);
		hardness = 0.5 * (hardness_F1 + hardness_F2);

		if (F1 == 6) {
			// Large surface hoar stays longer unstable! 1 dec 2007 (sb)
			const double grain_size = 2.*Edata.rg;
			if (grain_size >= 5.) {
				hardness = 1.;
			} else {
				hardness = MIN(hardness, 2.);
			}
		}
	} else if (Edata.theta[ICE] <= 0.7) { // Melt-freeze crust MFcr
		const double res_water_cont = Edata.res_wat_cont;
		if (Edata.theta[WATER] < 0.3 * res_water_cont) {
			hardness = 5.;
		} else if (Edata.theta[WATER] < 0.6 * res_water_cont) {
			hardness = 4.5;
		} else if (Edata.theta[WATER] < 0.85 * res_water_cont) {
			hardness = 4.;
		} else {
			hardness = 3.;
		}
	} else { // Ice Formations IF
		hardness = 6.;
	}
	// Limit to range {1, 6}
	hardness = MAX(1., MIN(6., hardness));
	return hardness;
}

/**
 * @brief Assign hand hardness to snow types according to density and grain size, original Canadian version
 * @author Implemented by C. Fierz: Regression from ASARC database by Bruce Jamieson on 2002-08-14
 * @param Edata
 * @return hand hardness index (1)
 */
double Stability::setHandHardnessASARC(const ElementData& Edata)
{
	const double gsz = 2.*Edata.rg;
	double A=0., B=0., C=0.;

	int F1, F2, F3;
	typeToCode(&F1, &F2, &F3, Edata.type); // Decompose type in its constituents

	// all types except MFcr
	if( Edata.mk%100 < 20 ) {
		switch ( F1 ) {
			case 0: { // Graupel PPgp; empirical!
				A = 1.5;
				B = 0.;
				C = 0.;
				break;
			}
			case 1: { // Precipitation Particles PP
				A = 0.45;
				B = 0.0068;
				C =  0.;
				break;
			}
			case 2: { // Decomposing and Fragmented precipitation particles DF
				A =  0.;
				B = 0.0140;
				C =  0.;
				break;
			}
			case 3: { // Rounded Grains RG
				A =  1.94;
				B = 0.0073;
				C = -0.192;
				break;
			}
			case 4: {
				if ( F2 != 9 ) { // Faceted Crystals FC
					A =  0.;
					B = 0.0138;
					C = -0.284;
				} else { // Rounding faceted particles FCxr, because F1=9 does not occur in SNOWPACK
					A =  1.29;
					B = 0.0094;
					C = -0.350;
				};
				break;
			}
			case 5: {
				if ( gsz > 1.5 ) { // Depth hoar DH, small dataset (N=41) !!!
					A = -0.80;
					B = 0.0150;
					C = -0.140;
				} else { // use FC values for small depth hoar grains
					A =  0.00;
					B = 0.0138;
					C = -0.284;
				};
				break;
			}
			case 6: { // Surface hoar SH; empirical: index 1 to 2 from HOAR_DENSITY_BURIED to 250 kg m-3
				A = 1. - hoar_density_buried/(250. - hoar_density_buried);
				B = 1./(250. - hoar_density_buried);
				C = 0.;
				break;
			}
			case 7: { // Melt Forms MF
				if ( Edata.theta[WATER] < 0.001 ) { // MF, dry
					A = 2.14;
					B = 0.0048;
					C =  0.;
				} else { // MF, wet (LWC > 0.)
					A = 3.00;
					B = 0.0000;
					C =  0.;
				};
				break;
			}
			case 8: { // Ice layer IFil
				A =  6.;
				B = 0.;
				C =  0.;
				break;
			}
			case 9: { // Rounding faceted particles FCxr
				A =  1.29;
				B = 0.0094;
				C = -0.350;
				break;
			}
			default: {
				A = Constants::undefined;
				B = 0.;
				C = 0.;
				break;
			}
		}
	} else if (Edata.theta[ICE] <= 0.7) { // Melt-freeze crust MFcr
		if (Edata.theta[WATER] < 0.3 * Edata.res_wat_cont) {
			A = 5.;
		} else if (Edata.theta[WATER] < 0.6 * Edata.res_wat_cont) {
			A = 4.5;
		} else if (Edata.theta[WATER] < 0.85 * Edata.res_wat_cont) {
			A = 4.;
		} else {
			A = 3.;
		}
	} else { // Ice Formations IF
		A = 6.;
	}

	double hardness = A + B*Edata.Rho + C*gsz;
	if (F1 == 6) {
		hardness = MIN(hardness, 2.);
	}
	// Limit to range {1, 6}
	hardness = MAX(1., MIN(6., hardness));
	return(hardness);
}

/*
 * START OF STABILITY SECTION
*/
/**
 * @brief Returns the critical stress state of a layer given the temperature and plastic strain rate.
 * @param epsNeckDot Neck strain rate (s-1)
 * @param Ts Temperature of layer (K)
 * @return Critical stress (Pa)
 */
double Stability::compCriticalStress(const double& epsNeckDot, const double& Ts)
{
	const double sigBrittle=1.e7;   // Brittle fracture stress of ice (Pa)
	const double C1=-6.6249;     // Constant
	const double C2=6.0780e-2;   // Constant
	const double C3=-1.3380e-4;  // Constant
	const double P1=70.000;      // Constant (Pa)

	// Find the rate dependent friction angle phi
	const double epsa = fabs(epsNeckDot); // Absolute value of plastic strain rate
	const double phi = P1*pow(epsa, 0.23)*mio::Cst::to_rad; // Function of strain rate dependent failure surface

	// Hydrostatic melting pressure
	// NOTE this function returns negative values for
	//   Ts <=181.2 K and Ts >= 273.15 K.
	//   The argument to the square root below becomes
	//   negative for Ts <= 180.4 K and Ts >= 274 K
	//   The maximum of the function is reached at 227.2 K
	//   HACK use this value for temperatures below 227.2 K (Quick and dirty fix;-)
	const double temp = (Ts >= 227.2)? Ts : 227.2;
	const double Pm = (C1 + C2*temp + C3*Optim::pow2(temp)) * 1.e9;

	// Return the critical stress. TODO check that argument of sqrt is correctly written
	return (Pm * tan(phi) * sqrt(1. - (Pm/(Pm + sigBrittle))));
}

/**
 * @brief Returns the layer stability index
 * The intra-layer stability criteria is given by the ratio S_f = S_c/S_n where
 * S_n is the neck stress and S_c is the critical stress.  The critical stress is determined
 * in the function st_CriticalStress. This function might get a little more involved as
 * time goes on.
 * @param *Edata
 * @return Deformation rate index
 */
double Stability::setDeformationRateIndex(ElementData& Edata)
{
	// If you have less than 5% ice then say you know you have something unstable
	if ( Edata.theta[ICE] < 0.05 ) {
		return(0.1);
	}

	const double eps1Dot = 1.76e-7; // Unit strain rate (at stress = 1 MPa) (s-1)
	const double sig1 = 0.5e6;      // Unit stress from Sinha's formulation (Pa)
	const double sig = -Edata.C;   // Overburden stress, that is, absolute value of Cauchy stress (Pa)
	const double Te = MIN(Edata.Te, Edata.melting_tk); // Element temperature (K)

	// First find the absolute neck stress
	const double sigNeck = Edata.neckStressEnhancement() * (sig); // Neck stress (Pa)
	// Now find the strain rate in the neck
	const double epsNeckDot =  eps1Dot * SnLaws::snowViscosityTemperatureTerm(Te) * mio::Optim::pow3(sigNeck/sig1); // Total strain rate in the neck (s-1) NOTE is it used here only?
	// Return the stability index
	return (MAX(0.1, MIN(compCriticalStress(epsNeckDot, Te) / sigNeck, 6.)));
}

/**
 * @brief Initializes stability parameters
 * @param STpar
 * @param Xdata
 * @param SIdata
 * @param i_psi_ref Reference slope angle (deg)
 */
void Stability::initStability(const double& i_psi_ref, StabilityData& STpar,
                              SnowStation& Xdata, std::vector<InstabilityData>& SIdata)
{
	const size_t nN = Xdata.getNumberOfNodes();

	STpar.Sig_c2 = Constants::undefined;
	STpar.strength_upper = 1001.;
	STpar.psi_ref = i_psi_ref*mio::Cst::to_rad;
	STpar.cos_psi_ref = cos(STpar.psi_ref);
	STpar.sin_psi_ref = sin(STpar.psi_ref);
	STpar.sig_n = Constants::undefined;
	STpar.sig_s = Constants::undefined;
	STpar.alpha_max_rad = 54.3*mio::Cst::to_rad; // alpha_max(38.) = 54.3 deg (J. Schweizer, IB 712, SLF)

	for(size_t n=Xdata.SoilNode; n<nN; n++) {
		SIdata[n].ssi      = Stability::max_stability;
		Xdata.Ndata[n].S_n = Stability::max_stability;
		Xdata.Ndata[n].S_s = Stability::max_stability;
		if (n < nN-1) {
			Xdata.Edata[n].S_dr = Stability::max_stability;
		}
	}

	Xdata.S_d = Xdata.S_n = Xdata.S_s = Xdata.S_4 = Xdata.S_5 = Stability::max_stability;
	Xdata.z_S_d = Xdata.z_S_n = Xdata.z_S_s = Xdata.z_S_4 = Xdata.z_S_5 = 0.;
	Xdata.S_class1 = Xdata.S_class2 = -1;
}

/**
 * @brief Returns the skier's penetration depth Pk
 * Adapted from Jamieson & Johnston, Ann. Glaciol., 26, 296-302 (1998)
 * @param *Xdata
 */
double Stability::compPenetrationDepth(const SnowStation& Xdata)
{
	double rho_Pk = Constants::eps2, dz_Pk = Constants::eps2; // Penetration depth Pk, from mean slab density
	double top_crust = 0., thick_crust = 0.;  // Crust properties
	bool crust = false;                       // Checks for crust
	size_t e_crust = Constants::stundefined;

	const double cos_sl = Xdata.cos_sl; // Cosine of slope angle
	size_t e = Xdata.getNumberOfElements(); //HACK is this right? It should be nNodes+1
	while ((e-- > Xdata.SoilNode) && ((Xdata.cH - (Xdata.Ndata[e].z + Xdata.Ndata[e].u))/cos_sl < 0.3)) {
		rho_Pk += Xdata.Edata[e].Rho*Xdata.Edata[e].L;
		dz_Pk  += Xdata.Edata[e].L;
		// Test for strong mf-crusts MFcr.
		// Look for the first (from top) with thickness perp to slope > 3cm
		if (!crust) {
			if ( (Xdata.Edata[e].mk%100 >= 20) && (Xdata.Edata[e].Rho > 500.) ) {
				if (e_crust == Constants::stundefined) {
					e_crust = e;
					top_crust = (Xdata.Ndata[e+1].z + Xdata.Ndata[e+1].u)/cos_sl;
					thick_crust += Xdata.Edata[e].L;
				} else if ( ((e_crust - e) < 2) ) {
					thick_crust += Xdata.Edata[e].L;
					e_crust = e;
				}
			} else if (e_crust > 0) {
				if (thick_crust > Stability::min_thick_crust) {
					crust = true;
				} else {
					e_crust = Constants::stundefined;
					top_crust = 0.;
					thick_crust = 0.;
				}
			}
		}
	}
	rho_Pk /= dz_Pk; //weighted average density of the snow slab penetrated by the skier

	// NOTE Pre-factor 0.8 introduced May 2006 by S. Bellaire
	return MIN(0.8 * 43.3 / rho_Pk, ((Xdata.cH / cos_sl) - top_crust));
}

/**
 * @brief Computes normal and shear stresses (kPa) reduced to psi_ref
 * @param STpar
 * @param stress Overload perpendicular to slope (Pa)
 * @param cos_sl Cosine of slope angle (1)
 */
void Stability::compReducedStresses(const double& stress, const double& cos_sl, StabilityData& STpar)
{
	STpar.sig_n = -stress*mio::Optim::pow2(STpar.cos_psi_ref/cos_sl)/1000.;
	STpar.sig_s = (STpar.sig_n)*STpar.sin_psi_ref/STpar.cos_psi_ref;
}

/**
 * @brief DEFAULT: Estimates the critical shear stress based on appropriate parameterisations
 * @param *Edata Xdata->Edata[e+1]
 * @param *Ndata Xdata->Ndata[e+1]
 * @param STpar
 * @param cH Computed height of snow (m)
 * @param cos_sl Cosine of slope angle (1)
 * @param date
 * @return return false on error, true otherwise
 */
bool Stability::setShearStrengthDEFAULT(const double& cH, const double& cos_sl, const mio::Date& date,
                                        ElementData& Edata, NodeData& Ndata, StabilityData& STpar)
{
	bool prn_wrn = false; //turn to true to print warnings

	const double rho_ri = Edata.Rho/Constants::density_ice; // Snow density relative to ice
	int F1, F2, F3; // Grain shape
	typeToCode(&F1, &F2, &F3, Edata.type); // Determine majority grain shape

	// Determine critical shear stress of element (kPa)
	// 1. Conway
	const double Sig_cC = 19.5*rho_ri*rho_ri;

	// 2. Grain Type dependent mostly from Jamieson & Johnston,
	//    Ann. Glaciol., 26, 296-302 (1998) and Ann. Glaciol., 32, 59-69 (2001)
	double phi = 0.; // Normal load correction
	double Sig_c2 = -1.0; // Critical shear stress (kPa)
	double Sig_c3 = -1.0; // Critical shear stress (kPa)

	switch( F1 ) {
		case 0: // Graupel, from O. Abe, Ann. Glaciol., 38, (2004), size-effect corrected
			Sig_c2 = 0.65*(82.*pow(rho_ri, 2.8));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 1: // PP
			Sig_c2 = 2.85*exp(1.13*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.056 + 0.022*STpar.sig_n;
			break;
		case 2: // DF
			Sig_c2 = 8.75*exp(1.54*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 3: // RG
			Sig_c2 = 7.39*exp(1.20*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 6: // SH
			switch(Stability::sh_mod) {
				case 0: // original T. Chalmers
					Sig_c2 = 0.336 + 0.0139*(date.getJulian() - Edata.depositionDate.getJulian()) +
							1.18*STpar.sig_n/Optim::pow2(STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*IOUtils::K_TO_C(Edata.Te) + 0.0204*Edata.rg;
					break;
				case 1: // original T. Chalmers & accounting for Emin as 2*rg (ml 13 Feb 2003)
					Sig_c2 = 0.336 + 0.0139*(date.getJulian() - Edata.depositionDate.getJulian()) +
							1.18*STpar.sig_n/Optim::pow2(STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*IOUtils::K_TO_C(Edata.Te) + 0.0204*2.*Edata.rg;
					break;
				case 2: // New regression by Bruce Jamieson w/o Emin (14 Feb 2003)
					Sig_c2 = 0.429 + 0.0138*(date.getJulian() - Edata.depositionDate.getJulian()) +
							1.12*STpar.sig_n/Optim::pow2(STpar.cos_psi_ref) - 0.596*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0785 *
							cH/cos_sl - 27.1*Edata.L/cos_sl +
							0.0202*IOUtils::K_TO_C(Edata.Te);
					break;
				default:
					Sig_c2 = 1.0;
					break;
			}
			Sig_c2 = MAX(0.1, Sig_c2);
			Sig_c3 = 84.*exp(2.55*log(rho_ri));
			break;
		case 7: // MF
			Sig_c2 = 21.*exp(1.24*log(rho_ri));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		default: // FC, DH, FCmx
			Sig_c2 = 18.5*exp(2.11*log(rho_ri));
			assert(STpar.sig_n>0.); //in a few cases, we have received sig_n<0 from the caller
			if (STpar.sig_n>0.)
				Sig_c3 = 1.36*exp(0.55*log(STpar.sig_n/STpar.cos_psi_ref));
			// phi = 0.08*Sig_c2 + 0.224;
			// Above correction not used by Jamieson & Johnston (1998), but considered by Lehning et al., Ann. Glaciol., 38, 331-338 (2004)
			break;
	}

		// Hack for MFCs
	if ( Edata.mk % 100 >= 20 ) {
		Sig_c2 = Sig_c3 = 4.;
	}

	// Final assignements
	STpar.Sig_c2 = MIN(Sig_c2, STpar.strength_upper);
	Edata.s_strength = Sig_c2;
	STpar.strength_upper = Sig_c2;
	STpar.phi = phi;

	// Warning message may be enabled for large differences in snow shear stength models
	if (prn_wrn
		    && (((fabs(Sig_c2-Sig_cC)/Sig_cC) > 10.) || ((Sig_c3 > 0.)
		        && ((fabs(Sig_c3-Sig_cC)/Sig_cC > 10.)))) ) {
		prn_msg( __FILE__, __LINE__, "wrn", date,"Large difference in Snow Shear Stength (type=%d)", F1);
		prn_msg(__FILE__, __LINE__, "msg-", Date(), "Conway: %lf Sig_c2: %lf Sig_c3: %lf\n", Sig_cC, Sig_c2, Sig_c3);
		return false;
	} else {
		return true;
	}
} // End setShearStrengthDEFAULT

/**
 * @brief STRENGTH_NIED: Estimates the critical shear stress based on appropriate parameterisations adapted for Japan
 * @param *Edata Xdata->Edata[e+1]
 * @param *Ndata Xdata->Ndata[e+1]
 * @param STpar
 * @param cH Computed height of snow (m)
 * @param cos_sl Cosine of slope angle (1)
 * @param date
 * @return return false on error, true otherwise
 */
bool Stability::setShearStrengthSTRENGTH_NIED(const double& cH, const double& cos_sl, const mio::Date& date,
                                              ElementData& Edata, NodeData& Ndata, StabilityData& STpar)
{
	int    F1, F2, F3;             // Grain shape
	double Sig_cC, Sig_c2, Sig_c3; // Critical shear stress (kPa)
	double Sig_ET, Sig_DH;         //NIED (H. Hirashima)
	double phi;                    // Normal load correction
	double rho_ri;                 // Snow density relative to ice
	bool prn_wrn = false;

	// Snow density relative to ice
	rho_ri = Edata.Rho/Constants::density_ice;
	// Determine majority grain shape
	typeToCode(&F1, &F2, &F3, Edata.type);

	// Determine critical shear stress of element (kPa)
	// 1. Conway
	Sig_cC = 19.5*rho_ri*rho_ri;
	// 2. Grain Type dependent mostly from Jamieson,
	//    Ann. Glaciol., 26, 296-302 (2001) and Ann. Glaciol., 32, 59-69 (1998)
	phi = 0.;
	Sig_c2 = -1.0;
	Sig_c3 = -1.0;
	switch ( F1 ) {
		case 0: // Graupel, from O. Abe, Ann. Glaciol. 38 (2004), size-effect corrected
			Sig_c2 = 0.65*(82.*pow(rho_ri, 2.8));
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 1: // PP //NIED (H. Hirashima)
			Sig_c2= 9.4*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
			phi = 0.08*Sig_c2 + 0.056 + 0.022*STpar.sig_n;
			break;
		case 2: case 3: // DF & RG //NIED (H. Hirashima)
			Sig_c2= 9.4*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
			phi = 0.08*Sig_c2 + 0.224;
			break;
		case 6: // SH
			switch(Stability::sh_mod) {
				case 0: // original T. Chalmers
					Sig_c2 = 0.336 + 0.0139*(date.getJulian() - Edata.depositionDate.getJulian()) +
							1.18*STpar.sig_n/Optim::pow2(STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*IOUtils::K_TO_C(Edata.Te) + 0.0204*Edata.rg;
					break;
				case 1: // original T. Chalmers & accounting for Emin as 2*rg (ml 13 Feb 2003)
					Sig_c2 = 0.336 + 0.0139*(date.getJulian() - Edata.depositionDate.getJulian()) +
							1.18*STpar.sig_n/Optim::pow2(STpar.cos_psi_ref) - 0.625*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0804 *
							cH/cos_sl - 28.7*Edata.L/cos_sl +
							0.0187*IOUtils::K_TO_C(Edata.Te) + 0.0204*2.*Edata.rg;
					break;
				case 2: // New regression by Bruce Jamieson w/o Emin (14 Feb 2003)
					Sig_c2 = 0.429 + 0.0138*(date.getJulian() - Edata.depositionDate.getJulian()) +
							1.12*STpar.sig_n/Optim::pow2(STpar.cos_psi_ref) - 0.596*(cH -
							(Ndata.z + Ndata.u))/cos_sl + 0.0785 *
							cH/cos_sl - 27.1*Edata.L/cos_sl +
							0.0202*IOUtils::K_TO_C(Edata.Te);
					break;
				default:
					Sig_c2 = 1.0;
					break;
			}
			Sig_c2 = MAX (0.1, Sig_c2);
			Sig_c3 = 84.*exp(2.55*log(rho_ri));
			break;
		case 7: // MF //NIED (H. Hirashima)
			Sig_c2= 4.97*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
			phi = 0.08*Sig_c2 + 0.224;
			break;
		default: // FC, DH, FCmx
			Sig_c2 = 18.5*exp(2.11*log(rho_ri));
			//Sig_c2 = 0.0391*exp(0.0141*Edata.theta[ICE]*Constants::density_ice); //NIED (H. Hirashima)
			Sig_c3 = 1.36*exp(0.55*log(STpar.sig_n/STpar.cos_psi_ref));
			// phi = 0.08*Sig_c2 + 0.224;
			// Above correction not used by Jamieson & Johnston (1998), but considered by Lehning et al., Ann. Glaciol., 38, 331-338 (2004)
			break;
	}

	// Hack for MFCs; not used by //NIED (H. Hirashima)

	// Final assignements
	//NIED (H. Hirashima)
	Sig_ET = 9.4*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.91)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
	Sig_DH = 2.3*0.0001*pow(Edata.theta[ICE]*Constants::density_ice,2.78)*exp(-0.235*Edata.theta[WATER]*100.)/1000.;
	Ndata.Sigdhf = Sig_ET - Edata.dhf*(Sig_ET - Sig_DH);
	Ndata.S_dhf = (Ndata.Sigdhf + phi*STpar.sig_n)/STpar.sig_s;
	// original SNOWPACK
	STpar.Sig_c2 = MIN(Sig_c2, STpar.strength_upper);
	Edata.s_strength = Sig_c2;
	STpar.strength_upper = Sig_c2;
	STpar.phi = phi;

	// Warning message may be enabled to warn for large differences in snow shear stength models
	if (prn_wrn
	        && (((fabs(Sig_c2-Sig_cC)/Sig_cC) > 10.) || ((Sig_c3 > 0.)
	            && ((fabs(Sig_c3-Sig_cC)/Sig_cC > 10.)))) ) {
		prn_msg( __FILE__, __LINE__, "wrn", date,"Large difference in Snow Shear Stength (type=%d)", F1);
		prn_msg(__FILE__, __LINE__, "msg-", Date(), "Conway: %lf Sig_c2: %lf Sig_c3: %lf\n", Sig_cC, Sig_c2, Sig_c3);
		return false;
	} else {
		return true;
	}
} // End setShearStrengthSTRENGTH_NIED

/**
 * @brief Returns the natural stability index Sn
 * The classic natural stability index Sn, that is, the ratio of shear stress to shear strength (static)
 * @param STpar
 */
double Stability::setNaturalStabilityIndex(const StabilityData& STpar)
{
	// Limit natural stability index to range {0.05, Stability::max_stability}
	return(MAX(0.05, MIN(((STpar.Sig_c2 + STpar.phi*STpar.sig_n)/STpar.sig_s), Stability::max_stability)));
}

/**
 * @brief Returns the skier stability index Sk reduced to psi_ref (usually 38 deg => Sk_38)
 * The classic skier stability index Sk(psi_ref), using P. Foehn's formula
 * (IAHS No162, 1987, p201) for the skier (load of 85 kg on 1.7 m long skis) induced shear stress.
 * @param depth_lay Depth of layer to investigate (m)
 * @param STpar
 */
double Stability::setSkierStabilityIndex(const double& depth_lay, const StabilityData& STpar)
{
	if ( depth_lay > Constants::eps ) {
		const double Alpha_max = STpar.alpha_max_rad;
		// Skier contribution to shear stress at psi_ref (in rad, corresponds usually to 38 deg)
		// about 0.1523 kPa / depth_lay at psi_ref = 38 deg and Alpha_max = 54.3 deg
		// double delta_sig = 2. * 0.5 * cos(Alpha_max) * Optim::pow2( sin(Alpha_max) ) * sin(Alpha_max + STpar.psi_ref);
		double delta_sig = 2. * (85.*Constants::g/1.7) * cos(Alpha_max) * Optim::pow2( sin(Alpha_max) ) * sin(Alpha_max + STpar.psi_ref);
		delta_sig /= Constants::pi *  depth_lay * STpar.cos_psi_ref; // in Pa
		delta_sig /= 1000.; // convert to kPa
		// Limit skier stability index to range {0.05, Stability::max_stability}
		return(MAX(0.05, MIN(((STpar.Sig_c2 + STpar.phi*STpar.sig_n)/(STpar.sig_s + delta_sig)), Stability::max_stability)));
	} else {
		return(Stability::max_stability); // strictly speaking, Sk is not defined
	}
}

/**
 * @brief Returns the structural stability index SSI
 * Adds one lemon to Sk for each structural instability found, presently hardness and grain size differences
 * above a given threshold
 * @param Edata_lower Xdata->Edata[e]
 * @param Edata_upper Xdata->Edata[e+1]
 * @param Sk Skier stability index Sk (Xdata->Ndata[e+1].S_s)
 * @param SIdata [e+1]
 * @return SIdata.ssi [e+1]
 */
double Stability::setStructuralStabilityIndex(const ElementData& Edata_lower, const ElementData& Edata_upper,
                                              const double& Sk, InstabilityData& SIdata)
{
	const double thresh_dhard=1.5, thresh_dgsz=0.5; // Thresholds for structural instabilities
	//const int nmax_lemon = 2; //Maximum number of structural instabilities looked at ("lemons")

	SIdata.n_lemon = 0;
	SIdata.dhard = fabs(Edata_lower.hard - Edata_upper.hard);
	if ( SIdata.dhard > thresh_dhard ) {
		SIdata.n_lemon++;
	}
	SIdata.dgsz = 2.*fabs(Edata_lower.rg - Edata_upper.rg);
	//double ref_gs= MIN (Edata_lower.rg,Edata_upper.rg);
	//SIdata.dgsz = (fabs(Edata_lower.rg - Edata_upper.rg))/(ref_gs);
	if ( SIdata.dgsz > thresh_dgsz ) {
		SIdata.n_lemon++;
	}
	// Skier Stability Index (SSI)
	SIdata.ssi = static_cast<double>(Stability::nmax_lemon - SIdata.n_lemon) + Sk;
	// Limit stability index to range {0.05, Stability::max_stability}
	SIdata.ssi = MAX(0.05, MIN (SIdata.ssi, Stability::max_stability));

	return SIdata.ssi;
}

/**
 * @brief Returns the Profile Stability Classification (Schweizer-Wiesinger Method)
 * @param Xdata
 * @return false if error, true otherwise
 */
bool Stability::classifyProfileStability(SnowStation& Xdata)
{
	// Dereference the element pointer containing micro-structure data
	const size_t nE = Xdata.getNumberOfElements();
	vector<ElementData>& EMS = Xdata.Edata;
	vector<NodeData>& NDS = Xdata.Ndata;

	// Initialize
	const double cos_sl = Xdata.cos_sl;

	// Classify only for Snowpacks thicker than Stability::minimum_slab (vertically)
	if ( (NDS[nE].z+NDS[nE].u)/cos_sl < Stability::minimum_slab ) {
		Xdata.S_class2 = 5;
		return true;
	}

	// First, find mean, maximum and minimum hardness
	double mH = 0., maxH = 0., minH = 7.; // Mean Hardness and Threshold
	for (size_t e = Xdata.SoilNode; e < nE; e++) {
		maxH = MAX (maxH, EMS[e].hard);
		minH = MIN (minH, EMS[e].hard);
		mH += EMS[e].hard;
	}
	mH /= (double)nE; //mean hardness of profile
	const double thH = MIN (0.5*mH, 1.); //threshold for critical transitions

	// Now make the classification for all critical transitions and keep track ....
	double h_Slab = EMS[nE-1].L/cos_sl;
	unsigned int count=0;
	int S = 5;
	for (size_t e = nE-2; e > Xdata.SoilNode; e--) {
		h_Slab += EMS[e].L/cos_sl;
		const double delta_H = EMS[e+1].hard - EMS[e].hard;
		if ( fabs(delta_H) > thH ) {
			count++;
			const size_t e_weak = (delta_H < 0.)? e+1 : e;

			// Decompose grain type to determine majority shape F1
			int F1, F2, F3; // Grain shape
			typeToCode(&F1, &F2, &F3, EMS[e_weak].type);

			// Remember that S is initialized to 5!!!
			// First consider wet weak layer
			if (EMS[e_weak].theta[WATER] > 0.75 * EMS[e_weak].res_wat_cont) {
				if ( EMS[e_weak].mk % 100 < 20 ) {
					S = 1;
				}
			} else { // Then do some stuff for dry snow
				double mH_u = 0.;
				for (size_t ii = e_weak; ii < nE; ii++) {
					mH_u += EMS[e].hard;
				}
				mH_u /= (double)(nE - e_weak);
				if ( mH > 2. ) {
					if ( delta_H < 0. ) {
						// Proposal Fz; (see original in version 7.4
						if ( (mH_u < 2.5) || (h_Slab < 0.7) ) {
							if ( (mH_u > 2.) && (h_Slab > 0.5) ) {
								S = MIN(S,4);
							} else {
								S = MIN(S,3);
							}
						} else {
							if ( minH > 2. ) {
								S = MIN(S,4);
							} else {
								S = MIN(S,3);
							}
						}
					}
				} else if ( mH < 1.5 ) {
					if ( (EMS[e_weak].rg > 0.5) && ((F1 > 3) && (F1 < 7)) ) {
						if ( (EMS[e_weak].rg > 0.75) && ((mH_u < 1.5) && (maxH < 2.5)) ) {
							S = 1;
						} else {
							S = MIN (S, 2);
						}
					} else {
						S = MIN (S, 3);
					}
				} else {
					S = MIN (S, 3);
				}
			} // end dry snow
		} // if weak layer found; also end of loop over elements
	} // end for

	if ( count == 0 ) {
		if ( mH > 2. ) {
			S = MIN (S, 4);
		} else if ( mH < 1.5) {
			if ( maxH > 2.3 ) {
				S = MIN (S, 2);
			} else {
				S = 1;
			}
		} else {
			S = MIN (S, 3);
		}
	}

	Xdata.S_class2 = S;

	if ( (S > 5) || (S < 1) ) {
		return false;
	} else {
		return true;
	}
}  // End classifyProfileStability

/**
 * @brief "Pattern recognition" of 10 profile types according to Schweizer, J. and M. Luetschg (2001).
 * Schweizer, J. and M. Luetschg, <i>Characteristics of human-triggered avalanches</i>, 2001, Cold Reg. Sci. Technol. 33(2-3): 147-162.
 * Note that analysis is done on vertical snow height.
 * @param *Xdata
 * @return false on error, true otherwise
 */
bool Stability::recognizeProfileType(SnowStation& Xdata)
{
	const size_t n_window=5;                              // Window half-width in number of elements
	const double L_base_0=0.2;
	const double min_hard=19.472, slope_hard=150.;        // Constants to compute reduced hardness,
	                                                      // (N) and (N m-1), respectively

	// cos of slope angle to convert height and thickness to vertical values
	const double cos_sl = Xdata.cos_sl;
	// Vertical snow depth
	const double cH = (Xdata.cH - Xdata.Ground)/cos_sl;

	// Check for snow profile shallower than 1.5*L_base_0 m (not classifiable)
	if ( cH <= 1.5*L_base_0 ) {
		Xdata.S_class1 = -1;
		return true;
	}

	const size_t nE_s = Xdata.getNumberOfElements() - Xdata.SoilNode; //number of snow elements

	// Dereference element and node pointers
	ElementData *EMS = &Xdata.Edata[0];
	vector<NodeData>& NDS = Xdata.Ndata;

	//temporary vectors
	vector<double> z_el(nE_s, 0.0);                            // Vertical element heigth (m)
	vector<double> L_el(nE_s, 0.0);                            // Vertical element thickness (m)
	vector<double> hard(nE_s, 0.0);                            // Hardness in N
	vector<double> red_hard(nE_s, 0.0);                        // Reduced hardness in N
	vector<double> deltaN(nE_s, 0.0);                          // Difference in hardness between layers in N

	// Absolute and reduced hardness profiles (N)
	for(size_t idx = nE_s; idx --> 0; ) { //because it is decremented before executing anything
		const size_t ii = idx+Xdata.SoilNode; //true element index
		z_el[idx] = ( (NDS[ii].z + NDS[ii].u) + (NDS[ii+1].z + NDS[ii+1].u) ) * .5 / cos_sl;
		L_el[idx] = EMS[ii].L/cos_sl;
		hard[idx] = min_hard*pow(EMS[ii].hard, 2.3607);
		red_hard[idx] = hard[idx] - (min_hard + slope_hard*(cH - z_el[idx]));
		if ( (unsigned)idx == nE_s-1 ) {
			deltaN[idx] = fabs(red_hard[idx] - min_hard);
		} else {
			deltaN[idx] = fabs(red_hard[idx] - red_hard[idx+1]);
		}
	}

	// Check for base strength (bottom L_base_0 m of snow cover)
	// not considering basal melt-freeze crusts
	double L_base = L_base_0;
	double L_sum = 0.;
	double mean_hard=0.,mean_gsz = 0.; // Means
	const double thresh_hard = 19.472*pow(4., 2.3607); // Hardness threshold (N)
	bool mf_base = (hard[0] > thresh_hard);
	size_t e = 0; //element index

	while ( L_sum <= L_base ) {
		if ( mf_base && (hard[e] < thresh_hard) ) {
			mf_base = false;
			L_base -= L_sum;
			L_sum = 0.;
			mean_hard = mean_gsz = 0.;
		}
		L_sum += L_el[e];
		mean_hard += L_el[e]*hard[e];
		mean_gsz += L_el[e]*(2.*EMS[e+Xdata.SoilNode].rg);
		e++;
	}
	// Averages
	mean_hard /= L_sum;
	mean_gsz /= L_sum;

	// Weak or strong base?
	const bool weak_base = ((mean_hard <= 275.) && (mean_gsz > 0.9));

	// Seek extremes over profile depth
	// Initialise
	size_t e_min = MIN(e, nE_s - 1); //e is >=0
	size_t e_el = e_min;
	size_t e_max = MIN(e_el + n_window, nE_s - 1);
	// Extremes and extremes' absolute heights
	double sum_red_hard = 0.;
	double red_hard_max = -9999.;
	double red_hard_min = 9999.;
	double deltaN_max = -9999.;
	double z_red_hard_min = 9999.;
	double z_red_hard_max = 9999.;
	double z_deltaN_max = 9999.;

	// First evaluation
	L_sum = 0.;
	for (e = e_min; e <= e_max; e++) {
		L_sum += L_el[e];
		sum_red_hard += L_el[e]*red_hard[e];
	}

	// Use window width of 2*n_window+1 elements
	while ( e_el <= (nE_s-1) ) {
		if ( (e_el - e_min) > n_window ) {
			L_sum -= L_el[e_min];
			sum_red_hard -= L_el[e_min]*red_hard[e_min];
			e_min++;
		}
		if ( (e_max < (nE_s-1)) && ((e_max - e_el) < n_window) ) {
			e_max++;
			L_sum += L_el[e_max];
			sum_red_hard += L_el[e_max]*red_hard[e_max];
		}
		// Find extremes ...
		if ( sum_red_hard/L_sum > red_hard_max ) {
			red_hard_max = sum_red_hard/L_sum;
			z_red_hard_max = z_el[e_el];
		}
		if ( sum_red_hard/L_sum < red_hard_min ) {
			red_hard_min = sum_red_hard/L_sum;
			z_red_hard_min = z_el[e_el];
		}
		e_el++;
	}

	// Find extremes for deltaN (no window required)
	e = 0;
	while ( (z_el[e] < L_base_0) && (e < (nE_s-1)) ) {
		e++;
	}
	for (; e < (nE_s-1); e++) {
		if ( deltaN[e] > deltaN_max ) {
			deltaN_max = deltaN[e];
			z_deltaN_max = z_el[e];
		}
	}

	if ( !((red_hard_max > (-150.*cH)) && (red_hard_min < 1500.)) ) {
		Xdata.S_class1 = -1;
		return false;
	}

	// Classify
	// Max. Hardness at position of maximum
	const double hard_max = red_hard_max + (min_hard + slope_hard*(cH - z_red_hard_max));
	int prf_type=-1; // Profile type
	if ( weak_base ) {
		// Position of extremes
		const double pos_max = (z_red_hard_max - L_base)/(cH - L_base);
		// Assign weak profile type
		if ( red_hard_max < 50. ) {
			prf_type = 1;
		} else if ( pos_max <= 0.3 ) {
			prf_type = 4;
		} else if ( pos_max <= 0.7 ) {
			prf_type = 3;
		} else if ( pos_max <= 0.9 ) {
			prf_type = 2;
		} else if ( (pos_max) > 0.9 && (hard_max > thresh_hard) ) {
			prf_type = 5;
		} else {
			prf_type = 4;
		}
	} else {// strong base
		// Position of extremes
		const double pos_max = (z_red_hard_max - L_base)/(cH - L_base);
		const double pos_min = (z_red_hard_min - L_base)/(cH - L_base);
		const double pos_max_deltaN = (z_deltaN_max - L_base)/(cH - L_base);

		// Assign strong profile type
		if ( (pos_max_deltaN > 0.85) && (hard_max > thresh_hard) ) {
			prf_type = 9;
		} else if ( (deltaN_max > 150.) && (pos_max_deltaN > pos_min) ) {
			prf_type = 7;
		} else if ( pos_max < 0.3 ) {
			prf_type = 6;
		} else if ( hard_max > thresh_hard ) {
			if ( fabs(red_hard_max - red_hard_min) < 50. ) {
				prf_type = 10;
			} else if ( red_hard_max < 50. ) {
				prf_type = 8;
			} else {
				prf_type =6;
			}
		} else {
			prf_type = 0;
		}
	}
	// end of classify

	Xdata.S_class1 = prf_type;

	return true;
} // End of recognizeProfileType

/*
 *  END OF CLASSIFICATION SECTION
*/

/**
 * @brief On a beautiful morning in September, with Foehn winds outside and incredibly fresh
 * colors, Michael finally started to implement the Stability thing. He is not con-
 * vinced that it is going to work but still dreams of lying in the sun at the Davos
 * lake. The stability information will be based on a very empirical principle. First
 * a distinction is made between "direct action" and "slab" situations. The former
 * have to do with strain weakening during heavy snowfalls or during melt situations.
 * The original Bob intra-layer stability will be adapted for this situation and
 * complemented by the Conway approach. The latter will be handled by an adaptation
 * of the Schweizer - Wiesinger profile classification in combination with a more
 * conventional stability index based on critical shear strength values. Halleluja.
 * @param Mdata CurrentMeteo
 * @param Xdata Profile
 */
void Stability::checkStability(const CurrentMeteo& Mdata, SnowStation& Xdata)
{
	const double cos_sl = Xdata.cos_sl; // Cosine of slope angle

	// Dereference the element pointer containing micro-structure data
	const size_t nN = Xdata.getNumberOfNodes();
	const size_t nE = nN-1;
	vector<NodeData>& NDS = Xdata.Ndata;
	vector<ElementData>& EMS = Xdata.Edata;

	vector<InstabilityData> SIdata = vector<InstabilityData>(nN); // Parameters for structural instabilities
	StabilityData  STpar;        // Stability parameters

	initStability(Stability::psi_ref, STpar, Xdata, SIdata);
	if ( (nE < Xdata.SoilNode+1) || plastic ) { // Return if bare soil or PLASTIC
		return;
	}

	const double Pk = compPenetrationDepth(Xdata); // Skier penetration depth
	size_t e = nE; // Counter
	while (e-- > Xdata.SoilNode) {
		EMS[e].hard = CALL_MEMBER_FN(*this, mapHandHardness[hardness_parameterization])(EMS[e]);
		EMS[e].S_dr = setDeformationRateIndex(EMS[e]);
		compReducedStresses(EMS[e].C, cos_sl, STpar);

		if ( !(CALL_MEMBER_FN(*this, mapShearStrength[strength_model])(Xdata.cH, cos_sl, Mdata.date,
		                                                               EMS[e], NDS[e+1], STpar))) {
			prn_msg(__FILE__, __LINE__, "msg-", Date(), "Node %03d of %03d", e+1, nN);
		}
		NDS[e+1].S_n = setNaturalStabilityIndex(STpar);
		const double depth_lay = (Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl - Pk; // corrected for skier penetration depth Pk.
		NDS[e+1].S_s = setSkierStabilityIndex(depth_lay, STpar);
		if (e < nE-1)
			NDS[e+1].ssi = setStructuralStabilityIndex(EMS[e], EMS[e+1], NDS[e+1].S_s, SIdata[e+1]);
		else
			NDS[nN-1].ssi = SIdata[nN-1].ssi = Stability::max_stability;
	}

	// Now find the weakest point in the stability profiles for natural and skier indices
	// Initialize
	size_t Swl_lemon = 0; // Lemon counter
	double Swl_d, Swl_n, Swl_ssi, zwl_d, zwl_n, zwl_ssi; // Temporary weak layer markers
	double Swl_Sk38, zwl_Sk38;       // Temporary weak layer markers
	Swl_d = Swl_n = Swl_ssi = Swl_Sk38 = INIT_STABILITY;
	zwl_d = zwl_n = zwl_ssi = zwl_Sk38 = Xdata.cH;

	// Natural and "deformation rate" Stability Index
	// Discard Stability::minimum_slab (in m) at surface
	e = nE;
	while ((e-- > Xdata.SoilNode) && (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < Stability::minimum_slab)) {};
	if (e==static_cast<size_t>(-1)) e=0; //HACK: this is ugly: e got corrupted if SoilNode==0
	
	if ((e > Xdata.SoilNode) && (e != IOUtils::unodata)) {
		// Slab must be thicker than Stability::ground_rough (m)  for an avalanche to release.
		while ((e-- > Xdata.SoilNode) && ((NDS[e+1].z + NDS[e+1].u)/cos_sl > Stability::ground_rough)) {
			// "deformation rate" Stability Index: find minimum ...
			if (Swl_d > EMS[e].S_dr) {
				Swl_d = EMS[e].S_dr;
				zwl_d = (NDS[e].z + NDS[e+1].z + NDS[e].u + NDS[e+1].u)/2.;
			}
			// Natural Stability Index: find minimum ...
			if ( Swl_n > NDS[e+1].S_n ) {
				Swl_n = NDS[e+1].S_n;
				zwl_n = NDS[e+1].z + NDS[e+1].u;
			}
		}
		// Assign minimum to stability indices
		Xdata.S_d = Swl_d;    Xdata.z_S_d = zwl_d - Xdata.Ground;
		Xdata.S_n = Swl_n;    Xdata.z_S_n = zwl_n - Xdata.Ground;
	} else {
		// Assign bottom values to stability indices
		Xdata.S_d = EMS[Xdata.SoilNode].S_dr;  Xdata.z_S_d = EMS[Xdata.SoilNode].L;
		Xdata.S_n = NDS[Xdata.SoilNode+1].S_n; Xdata.z_S_n = EMS[Xdata.SoilNode].L;
	}

	// Skier Stability Index
	//   Snow depth must be larger than Stability::ground_rough (m) and at least Stability::min_depth_ssi (m)
	//   snow must be left after discarding Pk for a SSI value to be searched.
	if ((Xdata.cH/cos_sl > Stability::ground_rough) && ((Xdata.cH/cos_sl - Pk) > Stability::min_depth_ssi)) {
		// Discard penetration depth Pk (in m) at surface
		e = nE;
		while ((e-- > Xdata.SoilNode) && (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < Pk)) {};
		if (e==static_cast<size_t>(-1)) e=0; //HACK: this is ugly: e got corrupted if SoilNode==0
		
		if ((e > Xdata.SoilNode) && (e != IOUtils::unodata)) {
			// Only down to Pk + Stability::skier_depth (m)

			while ((e-- > Xdata.SoilNode) && (((Xdata.cH - (NDS[e+1].z + NDS[e+1].u))/cos_sl) < (Pk + Stability::skier_depth)) && ((NDS[e+1].z + NDS[e+1].u)/cos_sl > Stability::ground_rough)) {
				// Skier Stability Index: find minimum OR consider number of structural instabilities in case of near equalities

				if ( (Swl_ssi > SIdata[e+1].ssi) || ((fabs(Swl_ssi - SIdata[e+1].ssi) < 0.09) && (SIdata[e+1].n_lemon > Swl_lemon)) ) {
					Swl_ssi = SIdata[e+1].ssi;
					zwl_ssi = NDS[e+1].z + NDS[e+1].u ;
					Swl_lemon = SIdata[e+1].n_lemon;
					Swl_Sk38 = NDS[e+1].S_s;
					zwl_Sk38 = NDS[e+1].z + NDS[e+1].u;
				}
			}
			// Assign minimum to stability indices
			Xdata.S_s = Swl_Sk38; Xdata.z_S_s = zwl_Sk38 - Xdata.Ground;
			Xdata.S_4 = Swl_ssi;  Xdata.z_S_4 = zwl_ssi - Xdata.Ground;
		} else {
			// Assign bottom values to stability indices
			Xdata.S_s = NDS[Xdata.SoilNode+1].S_s; Xdata.z_S_s = EMS[Xdata.SoilNode].L;
			Xdata.S_4 = SIdata[Xdata.SoilNode+1].ssi; Xdata.z_S_4 = EMS[Xdata.SoilNode].L;
		}
	} else {
		// Assign top values to stability indices
		Xdata.S_s = Stability::max_stability; Xdata.z_S_s = Xdata.cH;
		Xdata.S_4 = SIdata[nN-1].ssi; Xdata.z_S_4 = Xdata.cH;
	}

	switch (Stability::prof_classi) {
		case 0:
			// Classify in poor, fair and good based on master thesis of S. Bellaire (September 2005)
			if ((Swl_ssi > 0.) && (Swl_ssi < 100.)) {
				if (Swl_ssi >= 1.55) {
					Xdata.S_class2 = 5;
				} else if (Swl_ssi >= 1.25) {
					Xdata.S_class2 = 3;
				} else {
					Xdata.S_class2 = 1;
				}
			} else {
				Xdata.S_class2 = -1;
			}
			break;
		case 1:
		// Classify in poor, fair and good based on re-analysis by Schweizer/Bellaire (paper CRST 46 (2006) 52-59)
			if ((Swl_ssi > 0.) && (Swl_ssi < 100.)) {
				if ( Swl_Sk38 >= 0.45 ) {
					Xdata.S_class2 = 5;
				} else if ( Swl_ssi >= 1.32 ) {
						Xdata.S_class2 = 3;
					} else {
						Xdata.S_class2 = 1;
					}
			} else {
				Xdata.S_class2 = -1;
			}
			break;
		case 2:
			// Classify in poor, fair and good based on re-analysis after recalibration of settling,
			// Nov 2007(see rev 250/251)
			if ((Swl_ssi > 0.) && (Swl_ssi < 100.)) {
				if ( Swl_lemon >= 2 ) {
					Xdata.S_class2 = 1;
				} else if (Swl_lemon == 1) {
					if (Swl_Sk38 < 0.48) {
						Xdata.S_class2 = 1;
					} else {
						if (Swl_Sk38 < 0.71) {
							Xdata.S_class2 = 3;
						} else {
							Xdata.S_class2 = 5;
						}
					}
				} else {
					Xdata.S_class2 = 3;
				}
			} else {
				Xdata.S_class2 = -1;
			}
			break;
		case 3:
			// Classify in 5 classes based on ideas from Schweizer & Wiesinger
			if (!classifyProfileStability(Xdata)) {
				prn_msg( __FILE__, __LINE__, "wrn", Mdata.date,
					    "Profile classification failed! (classifyProfileStability)");
			}
			break;
	}

	if (classify_profile) {
		// Profile type based on "pattern recognition"; N types out of 10
		// We assume that we don't need it in Alpine3D
		if (!recognizeProfileType(Xdata)) {
			prn_msg( __FILE__, __LINE__, "wrn", Mdata.date, "Profile not classifiable! (recognizeProfileType)");
		}
	}
} // End checkStability
