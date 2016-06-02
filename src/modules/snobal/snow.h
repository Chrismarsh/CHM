#ifndef SNOW_H
#define SNOW_H
 
#include "envphys.h"

/*
 * Snow modeling.
 */
 
/* ----------------------------------------------------------------------- */
 
/*
 *  Constants
 */

/*
 *  thermal emissivity of snow
 */
#define SNOW_EMISSIVITY		0.98
 
/* ------------------------------------------------------------------------ */
 
/*
 *  Macros
 */

/*
 *  thermal conductivity of snow (J/(m sec K))
 *    (after Yen, 1965, see Anderson, 1976, pg. 31)
 *
 *	rho = snow density (kg/m^3)
 */
#define KTS(rho)        CAL_TO_J(0.0077*((rho)/1000.0)*((rho)/1000.0))

/*
 *  melt (kg/m^2)
 *
 *	Q = available energy (J/m^2)
 */
#define MELT(Q)         ( (Q) / LH_FUS(FREEZE) )

/*
 *  snow - liquid water mass fraction.
 *    from:   Davis, et. al. (1985)
 *
 *	mw = mass of liquid water (kg)
 *	ms = total mass of snow (kg)
 */
#define SNOH2O_MF(mw,ms)        ( mw / ms )

/*
 *  snow - liquid water volume fraction
 *    from:   Davis, et. al. (1985)
 *
 *	mw   = mass of liquid water (kg)
 *	ms   = total mass of snow (kg)
 *	rhos = density of snow (kg/m^3)
*/
#define SNOH2O_VF(mw,ms,rhos)   ( SNOH2O_MF(mw,ms) * (rhos/RHO_W0) )

/*
 * snow porosity
 *    from:   Davis, et. al. (1985)
 *
 *	mw   = mass of liquid water (kg)
 *	ms   = total mass of snow (kg)
 *	rhos = density of snow (kg/m^3)
 */
#define SNO_PORO(mw,ms,rhos)    \
        ( (RHO_ICE - (rhos * (1.0 - SNOH2O_MF(mw,ms)))) / RHO_ICE )

/*
 *  snow water saturation
 *    from:   Davis, et. al. (1985)
 *
 *	mw   = mass of liquid water (kg)
 *	ms   = total mass of snow (kg)
 *	rhos = density of snow (kg/m^3)
 */
#define SNO_SAT(mw,ms,rhos)     \
	( (RHO_ICE * rhos * SNOH2O_MF(mw,ms)) / ((RHO_ICE - rhos) * RHO_W0) )

/*
 *  water retained by snow at given saturation (see SNO_SAT)
 *
 *	d    = total depth of snow (m)
 *	rhos = density of snow (kg/m^3)
 *	sat  = snow saturation (see SNO_SAT)
 */
#define H2O_LEFT(d,rhos,sat)    \
        ( (sat * d * RHO_W0 * (RHO_ICE - rhos)) / RHO_ICE )

/*
 *  'dry' snow density (without H2O) at a given total snow density
 *  (with H2O) and snow saturation
 *
 *	rhos = total density of snow (kg/m^3)
 *	sat  = snow saturation (see SNO_SAT)
 */
#define DRY_SNO_RHO(rhos,sat)    \
	( ((rhos) - (sat) * RHO_W0) / (1 - (sat) * RHO_W0 / RHO_ICE) )

/* ------------------------------------------------------------------------ */
 
/*
 * Library functions.
 */

extern double	g_snow(double rho1, double rho2, double ts1, double ts2,
		       double ds1, double ds2, double pa);
extern double	g_soil(double rho, double tsno, double tg, double ds,
		       double dg, double pa);
extern double	new_tsno(double spm, double t0, double ccon);

/* ------------------------------------------------------------------------ */

#endif  /* SNOW_H */
