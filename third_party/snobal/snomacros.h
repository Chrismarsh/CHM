/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once
namespace snobalMacros
{

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
#define SNOW_EMISSIVITY        0.98

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
#define SNOH2O_MF(mw, ms)        ( mw / ms )

/*
 *  snow - liquid water volume fraction
 *    from:   Davis, et. al. (1985)
 *
 *	mw   = mass of liquid water (kg)
 *	ms   = total mass of snow (kg)
 *	rhos = density of snow (kg/m^3)
*/
#define SNOH2O_VF(mw, ms, rhos)   ( SNOH2O_MF(mw,ms) * (rhos/RHO_W0) )

/*
 * snow porosity
 *    from:   Davis, et. al. (1985)
 *
 *	mw   = mass of liquid water (kg)
 *	ms   = total mass of snow (kg)
 *	rhos = density of snow (kg/m^3)
 */
#define SNO_PORO(mw, ms, rhos)    \
        ( (RHO_ICE - (rhos * (1.0 - SNOH2O_MF(mw,ms)))) / RHO_ICE )

/*
 *  snow water saturation
 *    from:   Davis, et. al. (1985)
 *
 *	mw   = mass of liquid water (kg)
 *	ms   = total mass of snow (kg)
 *	rhos = density of snow (kg/m^3)
 */
#define SNO_SAT(mw, ms, rhos)     \
    ( (RHO_ICE * rhos * SNOH2O_MF(mw,ms)) / ((RHO_ICE - rhos) * RHO_W0) )

/*
 *  water retained by snow at given saturation (see SNO_SAT)
 *
 *	d    = total depth of snow (m)
 *	rhos = density of snow (kg/m^3)
 *	sat  = snow saturation (see SNO_SAT)
 */
#define H2O_LEFT(d, rhos, sat)    \
        ( (sat * d * RHO_W0 * (RHO_ICE - rhos)) / RHO_ICE )

/*
 *  'dry' snow density (without H2O) at a given total snow density
 *  (with H2O) and snow saturation
 *
 *	rhos = total density of snow (kg/m^3)
 *	sat  = snow saturation (see SNO_SAT)
 */
#define DRY_SNO_RHO(rhos, sat)    \
    ( ((rhos) - (sat) * RHO_W0) / (1 - (sat) * RHO_W0 / RHO_ICE) )

/*
 *  Stefan-Boltzmann constant (W / m^2 / deg^4)
 */
#define STEF_BOLTZ      5.67032e-8

/*
 *  Planck radiation physics
 */
#define hPLANCK         6.626176e-34    /*  J sec                       */
#define PLANCK1st       3.741832e-16
#define PLANCK2nd       1.438786e-2
#define kBOLTZ          1.380662e-23    /*  J / deg                     */

/*
 *  Methods for computing gamma values (see the function 'mwgamma')
 */
#define DELTA_EDDINGTON        0
#define MEADOR_WEAVER        1

/* 	 
*  Maximum density due to compaction (kg/m^3). 	 
*/ 
#define	RMX	550 	

/* 	 
*  days over which compaction occurs 	 
*/ 
#define RD1	23.5 
#define RD2	24.5 

/* 	 
*  max swe to consider - after this swe value, compaction is maximized 	 
*/ 
#define SWE_MAX	2000.0 

/*
 * Status codes returned by integer functions.
 */


#define OK		1
#define ERROR		(-2)	/* from IPW version 1.0 */

#define AH        1.0    /* ratio sensible/momentum phi func	*/
#define AV        1.0    /* ratio latent/momentum phi func	*/
#define ITMAX      200  //50    /* max # iterations allowed		*/
#define PAESCHKE    7.35    /* Paeschke's const (eq. 5.3)		*/
//#define THRESH        1.e-1    /* convergence threshold	 1e-5 is orign	*/
const double THRESH = 1.e-1;

#define SM        0
#define SH        1
#define SV        2
#define BETA_S        5.2
#define BETA_U        16

/*
 * Environmental physics.
 */

/* ----------------------------------------------------------------------- */

/*
 *  Constants
 */
#define PI	3.14159265 

/*
 *  molecular weight of air (kg / kmole)
 */
#define MOL_AIR         28.9644

/*
 *  molecular weight of water vapor (kg / kmole)
 */
#define MOL_H2O         18.0153

/*
 *  gas constant (J / kmole / deg)
 */
#define RGAS            8.31432e3

/*
 *  triple point of water at standard pressure (deg K)
 */
#define FREEZE          2.7316e2
#define BOIL            3.7315e2

/*
 *  specific heat of air at constant pressure (J / kg / deg)
 */
#define CP_AIR          1.005e3

/*
 *  specific heat of water at 0C (J / (kg K))
 */
#define CP_W0    4217.7

/*
 *  density of water at 0C (kg/m^3)
 *    (from CRC handbook pg F-11)
 */
#define RHO_W0    999.87

/* 	
 *  Density of water (kg/m^3). 	
 *  Use by updated compaction parameterization.
 * TODO: coherence between water and RHO_W0
*/ 
#define water	1000.0 	


/*  density of ice - no air (kg/m^3)
 *    (from CRC handbook pg F-1)
 */
#define RHO_ICE    917.0

/*  thermal conductivity of wet sand (J/(m sec K))
 *    (from Oke, 1978, pg. 38)
 */

// set in the main class now

//#define KT_WETSAND      2.2

/*
 *  standard sea level pressure (Pa)
 */
#define SEA_LEVEL       1.013246e5

/*
 *  standard sea level air temp (K)
 */
#define STD_AIRTMP      2.88e2

/*
 *  standard lapse rate (K/m)
 */
#define STD_LAPSE_M      -0.0065

/*
 *  standard lapse rate (K/km)
 */
#define STD_LAPSE       -6.5

/*
 *  dry adiabatic lapse rate (deg / m)
 */
#define DALR        ( GRAVITY / CP_AIR )

/*
 *  gravitational acceleration at reference latitude 45d 32m 33s (m/s^2)
 */
#define GRAVITY        9.80665

/*
 *  Earth equivalent spherical radius (km)
 */
#define EARTH_RADIUS    6.37122e3

/*
 *  velocity of light, m/s (Marx, Nature, 296, 11, 1982)
 */
#define LIGHT_SPEED     2.99722458e8

/*
 *  Von Karman constant
 */
#define VON_KARMAN      3.5e-1

/* ------------------------------------------------------------------------ */

/*
 *  Macros
 */

/*
 *  equation of state, to give density of a gas (kg/m^3)
 *
 *	p = pressure (Pa)
 *	m = molecular weight (kg/kmole)
 *	t = temperature (K)
 *
 *  or, inversely, to give pressure (Pa)
 *
 *      rho = density (kg/m^3)
 *	m   = molecular weight (kg/kmole)
 *	t   = temperature (K)
 */
#define GAS_DEN(p, m, t)          ((p)*(m)/(RGAS*(t)))
#define EQ_STATE(rho, m, t)       ((rho)*(RGAS)*(t)/(m))

/*
 *  virtual temperature, i.e. the fictitious temperature that air must
 *  have at the given pressure to have the same density as a water vapor
 *  and air mixture at the same pressure with the given temperature and
 *  vapor pressure.
 *
 *      t = temperature (K)
 *      e = vapor pressure
 *	P = pressure (e and P in same units),
 */
#define VIR_TEMP(t, e, P)        ((t)/(1.-(1.-MOL_H2O/MOL_AIR)*((e)/(P))))

/*
 *  inverse of VIR_TEMP
 *
 *      tv = virtual temperature (K)
 *      e  = vapor pressure
 *	P  = pressure (e and P in same units),
 */
#define INV_VIR_TEMP(tv, e, P)          ((tv)*(1.-(1.-MOL_H2O/MOL_AIR)*((e)/(P))))

/*
 *  potential temperature
 *
 *      t = temperature (K)
 *	p = pressure (Pa)
 */
#define POT_TEMP(t, p)           ((t)*pow(1.e5/(p),RGAS/(MOL_AIR*CP_AIR)))

/*
 *  inverse of POT_TEMP
 *
 *	theta = potential temperature (K)
 *	p     = pressure (Pa)
 */
#define    INV_POT_TEMP(theta, p)    ((theta)/pow(1.e5/(p),RGAS/(MOL_AIR*CP_AIR)))

/*
 *  specific heat of ice (J/(kg K))
 *    (from CRC table D-159; most accurate from 0 to -10 C)
 *
 *	t = temperature (K)
 */
#define CP_ICE(t)       ( CAL_TO_J(0.024928 + (0.00176*(t))) / G_TO_KG(1) )

/*
 *  specific heat of water (J/(kg K))
 *    (from CRC table D-158; most accurate from 0 to +10 C)
 *    (incorrect at temperatures above 25 C)
 *
 *	t = temperature (K)
 */
#define CP_WATER(t)     ( CP_W0 - 2.55*((t)-FREEZE) )

/*
 *  integral of hydrostatic equation over layer with linear temperature
 *  variation
 *
 *	pb = base level pressure
 *	tb = base level temp (K)
 *	L  = lapse rate (deg/km)
 *	h  = layer thickness (km)
 *      g  = grav accel (m/s^2)
 *	m  = molec wt (kg/kmole)
 *
 *	(the factors 1.e-3 and 1.e3 are for units conversion)
 */
#define HYSTAT(pb, tb, L, h, g, m)           ((pb) * (((L)==0.) ?\
                exp(-(g)*(m)*(h)*1.e3/(RGAS*(tb))) :\
                pow((tb)/((tb)+(L)*(h)),(g)*(m)/(RGAS*(L)*1.e-3))))

/*
 *  inverse of integral of hydrostatic equation over layer with linear
 *  temperature variation
 *
 *      pb = base level pressure
 *	tb = base level temp (K)
 *      hb = base level geopotential altitude (km)
 *      p  = level pressure
 *	t  = level temperature (K)
 *      g  = grav accel (m/s^2)
 *	m  = molec wt (kg/kmole)
 *
 *      (the factor 1.e-3 is for units conversion)
 */
#define INV_HYSTAT(pb, tb, hb, p, t, g, m) ((hb)+\
                    1.e-3*log((p)/(pb))*(RGAS/((g)*(m)))*\
                    (((tb)==(t)) ? (-(t)) :\
                    ((t)-(tb))/log((tb)/(t))))

/*
 *  specific humidity
 *
 *	e = vapor pressure
 *	P = pressure (same units as e)
 */
#define SPEC_HUM(e, P)     ((e)*MOL_H2O/(MOL_AIR*(P)+(e)*(MOL_H2O-MOL_AIR)))

/*
 *  vapor pressure from specific humidity
 *
 *	q = specific humidity
 *	P = pressure
 */
#define INV_SPEC_HUM(q, P)  (-MOL_AIR*(P)*(q)/((MOL_H2O-MOL_AIR)*(q)-MOL_H2O))

/*
 *  mixing ratio
 *
 *	e = vapor pressure
 *	P = pressure (same units as e)
 */
#define MIX_RATIO(e, P)          ((MOL_H2O/MOL_AIR)*(e)/((P)-(e)))

/*
 *  vapor pressure from mixing ratio
 *
 *	w = mixing ratio
 *	P = pressure
 */
#define INV_MIX_RATIO(w, P)    ((w)*(P)/((w)+(MOL_H2O/MOL_AIR)))

/*
 *  latent heat of vaporization
 *
 *	t = temperature (K)
 */
#define LH_VAP(t)               (2.5e6 - 2.95573e3 *((t) - FREEZE))

/*
 *  latent heat of fusion
 *
 *	t = temperature (K)
 */
#define LH_FUS(t)               (3.336e5 + 1.6667e2 * (FREEZE - (t)))

/*
 *  latent heat of sublimination (J/kg)
 *    from the sum of latent heats of vaporization and fusion,
 *
 *	t = temperature (K)
 */
#define LH_SUB(t)        ( LH_VAP(t) + LH_FUS(t) )


/*
 *  effectuve diffusion coefficient (m^2/sec) for saturated porous layer
 *  (like snow...).  See Anderson, 1976, pg. 32, eq. 3.13.
 *
 *	pa = air pressure (Pa)
 *	ts = layer temperature (K)
 */
#define DIFFUS(pa, ts)   ( (0.65*(SEA_LEVEL/(pa)) * \
                        pow(((ts)/FREEZE),14.0)) * (0.01*0.01) )

/*
 *  water vapor flux (kg/(m^2 sec)) between two layers
 *
 *	air_d = air density (kg/m^3)
 *	k     = diffusion coef. (m^2/sec)
 *	q_dif = specific hum. diff between layers (kg/kg)
 *	z_dif = absolute distance between layers (m)
 *
 *	note:   q_dif controls the sign of the computed flux
 */
#define EVAP(air_d, k, q_dif, z_dif)       ( air_d * k * (q_dif/z_dif) )

/*
 *  dry static energy (J/kg)
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	t = air temperature (K)
 *	z = elevation (m)
 */
#define DSE(t, z)        ( (CP_AIR * (t) ) + (GRAVITY * (z) ) )

/*
 *  inverse of DSE: air temperature (K) from dry static energy
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	dse = dry static energy (J/kg)
 *	z   = elevation (m)
 */
#define INV_DSE(dse, z)  ((dse) - (GRAVITY * (z)) / CP_AIR )

/*
 *  moist static energy (J/kg)
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	z = elevation (m)
 *	t = air temperature (K)
 *	w = mixing ratio
 */
#define MSE(z, t, w)      ( DSE((t),(z)) + LH_VAP(t) * (w) )

/*
 *  inverse of MSE: vapor press (Pa) from moist static energy
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	z   = elevation (m)
 *	t   = air temperature (K)
 *	mse = moist static energy (J/kg)
 */
#define INV_MSE(z, t, mse)        ( ((mse) - DSE((t),(z))) / LH_VAP(t) )

/* ------------------------------------------------------------------------ */
#define MAX_SNOW_DENSITY	750

/*
*  default for snowcover's maximum liquid h2o content as volume
*  ratio: V_water/(V_snow - V_ice)
*/
#define DEFAULT_MAX_H2O_VOL  0.01

/*
 *  default for maximum active (surface) layer depth (m)
 */
#define DEFAULT_MAX_Z_S_0  0.25

/*
 *  default for depth of soil temperature measurement (m)
 */
#define DEFAULT_Z_G      0.5

/*
 *  Minimum valid snow temperature (C).  This is also what temperatures
 *  are set to when there's no snow (instead of 0 K).  This yields a
 *  smaller quantization range in the output image: -75 C to 0 C
 *  (instead of -273.16 C to 0 C).
 */
#define MIN_SNOW_TEMP    -75


/*
 *  default for medium run timestep (minutes)
 */
#define    DEFAULT_MEDIUM_TSTEP  15

/*
 *  default for small run timestep (minutes)
 */
#define    DEFAULT_SMALL_TSTEP  1

/*
 *  default for normal run timestep's threshold for a layer's mass
 *  (kg/m^2)
 */
#define    DEFAULT_NORMAL_THRESHOLD  60.0;

/*
 *  default for medium run timestep's threshold for a layer's mass
 *  (kg/m^2)
 */
#define    DEFAULT_MEDIUM_THRESHOLD 10.0

/*
 *  default for small run timestep's threshold for a layer's mass
 *  (kg/m^2)
 */
#define    DEFAULT_SMALL_THRESHOLD  1.0

/*
 *  Does a time fall within the current input data timestep?
 */
#define IN_CURR_DATA_TSTEP(time)    \
        ((current_time <= (time)) && \
         ((time) < current_time + tstep_info[DATA_TSTEP].time_step))

/*
 *  Units of physical measurement
 */

/* ------------------------------------------------------------------------ */

/*
 *  convert Celsius to Kelvin
 */
//#define C_TO_K(c)        ((c) + FREEZE)

/*
 *  convert Kelvin to Celsius
 */
//#define K_TO_C(k)        ((k) - FREEZE)

/*
 *  Convert kilograms to grams.
 */
#define KG_TO_G(kg)        ((kg) * 1000.0)

/*
 *  Convert grams to kilograms.
 */
#define G_TO_KG(g)        ((g) * 0.001)

/*
 *  Convert Joules to calories.
 */
#define J_TO_CAL(j)        ((j) * 0.238846)

/*
 *  Convert calories to Joules
 */
#define CAL_TO_J(c)        ((c) * 4.186798188)

/*
 *  convert wavelength (um) to wave number (1/cm)
 */
#define WAVENO(x)                       ( (int)(10000. / (x) + .5) )

/*
 *  convert wave number (1/cm) to wavelength (um)
 */
#define WAVELEN(nu)                     ( 10000. / (nu) )
/* 	 
*  seconds in an hour 	 
*/
#define nsec_hour	3600 	
}
