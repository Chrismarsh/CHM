#ifndef RADIATION_H
#define RADIATION_H
 
/*
 * Radiation modeling.
 */
 
/* ----------------------------------------------------------------------- */
 
/*
 *  Constants
 */

/*
 *  Stefan-Boltzmann constant (W / m^2 / deg^4)
 */
#define STEF_BOLTZ      5.67032e-8
 
/*
 *  Planck radiation constants
 */
#define hPLANCK         6.626176e-34    /*  J sec                       */
#define PLANCK1st       3.741832e-16
#define PLANCK2nd       1.438786e-2
#define kBOLTZ          1.380662e-23    /*  J / deg                     */

/*
 *  Methods for computing gamma values (see the function 'mwgamma')
 */
#define DELTA_EDDINGTON		0
#define MEADOR_WEAVER		1

/* ------------------------------------------------------------------------ */
 
/*
 *  Macros
 */

/* ------------------------------------------------------------------------ */
 
/*
 * Library functions.
 */

extern double	beta_0(double u0, double g);
extern double	brutsaert(double ta, double lambda, double ea, double z,
			  double pa);
extern void	delta_edd(double *omega, double *g, double *tau);
extern void	mwgamma(double u0, double w, double g, double *gam, int method);
extern double	net_therm(double ta, double ts, double pa, double ea,
			  double surfemiss, double skyfac, double lapse,
			  double z);
extern int	twostream(double *gamma, double omega, double mu0, double tau,
			  double r0, double *refl, double *trans,
			  double *btrans);


#endif  /* RADIATION_H */
