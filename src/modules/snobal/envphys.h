#ifndef ENVPHYS_H
#define ENVPHYS_H
 
#include <math.h>

#include "_physunits.h"

/*
 * Environmental physics.
 */
 
/* ----------------------------------------------------------------------- */
 
/*
 *  Constants
 */
 
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
#define CP_W0   	4217.7

/*
 *  density of water at 0C (kg/m^3)
 *    (from CRC handbook pg F-11)
 */
#define RHO_W0  	999.87

/*  density of ice - no air (kg/m^3)
 *    (from CRC handbook pg F-1)
 */
#define RHO_ICE 	917.0

/*  thermal conductivity of wet sand (J/(m sec K))
 *    (from Oke, 1978, pg. 38)
 */
#define KT_WETSAND      2.2

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
#define DALR		( GRAVITY / CP_AIR )

/*
 *  gravitational acceleration at reference latitude 45d 32m 33s (m/s^2)
 */
#define GRAVITY		9.80665

/*
 *  Earth equivalent spherical radius (km)
 */
#define EARTH_RADIUS	6.37122e3

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
#define GAS_DEN(p,m,t)          ((p)*(m)/(RGAS*(t)))
#define EQ_STATE(rho,m,t)       ((rho)*(RGAS)*(t)/(m))

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
#define VIR_TEMP(t,e,P)		((t)/(1.-(1.-MOL_H2O/MOL_AIR)*((e)/(P))))

/*
 *  inverse of VIR_TEMP
 *
 *      tv = virtual temperature (K)
 *      e  = vapor pressure
 *	P  = pressure (e and P in same units),
 */
#define INV_VIR_TEMP(tv,e,P)          ((tv)*(1.-(1.-MOL_H2O/MOL_AIR)*((e)/(P))))

/*
 *  potential temperature
 *
 *      t = temperature (K)
 *	p = pressure (Pa)
 */
#define POT_TEMP(t,p)           ((t)*pow(1.e5/(p),RGAS/(MOL_AIR*CP_AIR)))

/*
 *  inverse of POT_TEMP
 *
 *	theta = potential temperature (K)
 *	p     = pressure (Pa)
 */
#define	INV_POT_TEMP(theta,p)	((theta)/pow(1.e5/(p),RGAS/(MOL_AIR*CP_AIR)))

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
#define HYSTAT(pb,tb,L,h,g,m)           ((pb) * (((L)==0.) ?\
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
#define INV_HYSTAT(pb,tb,hb,p,t,g,m) ((hb)+\
	                1.e-3*log((p)/(pb))*(RGAS/((g)*(m)))*\
	       	        (((tb)==(t)) ? (-(t)) :\
	                ((t)-(tb))/log((tb)/(t))))

/*
 *  specific humidity
 *
 *	e = vapor pressure
 *	P = pressure (same units as e)
 */
#define SPEC_HUM(e,P)     ((e)*MOL_H2O/(MOL_AIR*(P)+(e)*(MOL_H2O-MOL_AIR)))

/*
 *  vapor pressure from specific humidity
 *
 *	q = specific humidity
 *	P = pressure
 */
#define INV_SPEC_HUM(q,P)  (-MOL_AIR*(P)*(q)/((MOL_H2O-MOL_AIR)*(q)-MOL_H2O))

/*
 *  mixing ratio
 *
 *	e = vapor pressure
 *	P = pressure (same units as e)
 */
#define MIX_RATIO(e,P)          ((MOL_H2O/MOL_AIR)*(e)/((P)-(e)))

/*
 *  vapor pressure from mixing ratio
 *
 *	w = mixing ratio
 *	P = pressure
 */
#define INV_MIX_RATIO(w,P)	((w)*(P)/((w)+(MOL_H2O/MOL_AIR)))

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
#define LH_SUB(t)       	( LH_VAP(t) + LH_FUS(t) )


/*
 *  effectuve diffusion coefficient (m^2/sec) for saturated porous layer
 *  (like snow...).  See Anderson, 1976, pg. 32, eq. 3.13.
 *
 *	pa = air pressure (Pa)
 *	ts = layer temperature (K)
 */
#define DIFFUS(pa,ts)   ( (0.65*(SEA_LEVEL/(pa)) * \
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
#define EVAP(air_d,k,q_dif,z_dif)       ( air_d * k * (q_dif/z_dif) )

/*
 *  dry static energy (J/kg)
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	t = air temperature (K)
 *	z = elevation (m)
 */
#define DSE(t,z)        ( (CP_AIR * (t) ) + (GRAVITY * (z) ) )

/*
 *  inverse of DSE: air temperature (K) from dry static energy
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	dse = dry static energy (J/kg)
 *	z   = elevation (m)
 */
#define INV_DSE(dse,z)  ((dse) - (GRAVITY * (z)) / CP_AIR )

/*
 *  moist static energy (J/kg)
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	z = elevation (m)
 *	t = air temperature (K)
 *	w = mixing ratio
 */
#define MSE(z,t,w)      ( DSE((t),(z)) + LH_VAP(t) * (w) )

/*
 *  inverse of MSE: vapor press (Pa) from moist static energy
 *    (pg. 332, "Dynamic Meteorology", Holton, J.R., 1979)
 *
 *	z   = elevation (m)
 *	t   = air temperature (K)
 *	mse = moist static energy (J/kg)
 */
#define INV_MSE(z,t,mse)        ( ((mse) - DSE((t),(z))) / LH_VAP(t) )

/* ------------------------------------------------------------------------ */
 
/*
 * Library functions.
 */

extern double	bevap(double netrad, double advec, double bowen,
		      double storage, double ts);
extern double	bowen(double p, double ta, double ts, double ea, double es);
extern int	budyer(double z, double z0, double t, double t0,
		       double e, double e0, double u, double u0, double  p,
		       double *h, double *le);
extern int	budyer2(double zu, double zt, double z0, double t, double t0,
			double e, double e0, double u, double p,
			double *h, double *le);
extern double	dew_point(double e);
extern double	efcon(double k, double t, double p);
extern double	evap(double le, double ts);
extern double	heat_stor(double cp, double spm, double tdif);
extern int      hle1(double press, double ta, double ts, double za,
		     double ea, double es, double zq, double u, double zu,
		     double z0, double *h, double *le, double *e);
extern double   psychrom(double tdry, double twet, double press);
extern double   ri_no(double z2, double z1, double t2, double t1,
		      double u2, double u1);
extern double	sati(double tk);
extern double	satw(double tk);
extern double	ssxfr(double  k1, double  k2, double  t1, double  t2,
		      double  d1, double  d2);

/* ------------------------------------------------------------------------ */

#endif  /* ENVPHYS_H */
