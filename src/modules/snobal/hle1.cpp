#include <math.h>

#include "ipw.h"
#include "envphys.h"

#define AH		1.0	/* ratio sensible/momentum phi func	*/
#define AV		1.0	/* ratio latent/momentum phi func	*/
#define ITMAX		50	/* max # iterations allowed		*/
#define PAESCHKE	7.35	/* Paeschke's const (eq. 5.3)		*/
#define THRESH		1.e-5	/* convergence threshold		*/

#define SM		0
#define SH		1
#define SV		2
#define BETA_S		5.2
#define BETA_U		16
#include <string>
void usrerr(std::string error, std::string error1){}
/* ----------------------------------------------------------------------- */

/*
 * psi-functions
 *	code =	SM	momentum
 *		SH	sensible heat flux
 *		SV	latent heat flux
 */

static double
psi(
	double	zeta,		/* z/lo				*/
	int	code)		/* which psi function? (see above) */
{
	double	x;		/* height function variable	*/
	double	result;

	if (zeta > 0) {		/* stable */
		if (zeta > 1)
			zeta = 1;
		result = -BETA_S * zeta;
	}

	else if (zeta < 0) {	/* unstable */

		x = sqrt(sqrt(1 - BETA_U * zeta));

		switch (code) {
		case SM:
			result = 2 * log((1+x)/2) + log((1+x*x)/2) -
				2 * atan(x) + M_PI_2;
			break;

		case SH:
		case SV:
			result = 2 * log((1+x*x)/2);
			break;

		default: /* shouldn't reach */
				break;
//			bug("psi-function code not of these: SM, SH, SV");
		}
	}

	else {			/* neutral */
		result = 0;
	}

	return (result);
}

/* ----------------------------------------------------------------------- */

int
hle1(
	double	press,	/* air pressure (Pa)			*/
	double	ta,	/* air temperature (K) at height za	*/
	double	ts,	/* surface temperature (K)		*/
	double	za,	/* height of air temp measurement (m)	*/
	double	ea,	/* vapor pressure (Pa) at height zq	*/
	double	es,	/* vapor pressure (Pa) at surface	*/
	double	zq,	/* height of spec hum measurement (m)	*/
	double	u,	/* wind speed (m/s) at height zu	*/
	double	zu,	/* height of wind speed measurement (m)	*/
	double	z0,	/* roughness length (m)			*/

    /* output variables */

	double *h,	/* sens heat flux (+ to surf) (W/m^2)	*/
	double *le,	/* latent heat flux (+ to surf) (W/m^2)	*/
	double *e)	/* mass flux (+ to surf) (kg/m^2/s)	*/
{
	double	ah = AH;
	double	av = AV;
	double	cp = CP_AIR;
	double	d0;	/* displacement height (eq. 5.3)	*/
	double	dens;	/* air density				*/
	double	diff;	/* difference between guesses		*/
	double	factor;
	double	g = GRAVITY;
	double	k = VON_KARMAN;
	double	last;	/* last guess at lo			*/
	double	lo;	/* Obukhov stability length (eq. 4.25)	*/
	double	ltsh;	/* log ((za-d0)/z0)			*/
	double	ltsm;	/* log ((zu-d0)/z0)			*/
	double	ltsv;	/* log ((zq-d0)/z0)			*/
	double	qa;	/* specific humidity at height zq	*/
	double	qs;	/* specific humidity at surface		*/
	double	ustar;	/* friction velocity (eq. 4.34')	*/
	double	xlh;	/* latent heat of vap/subl		*/
	int	ier;	/* return error code			*/
	int	iter;	/* iteration counter			*/

	/*
	 * check for bad input
	 */

	/* heights must be positive */
	if (z0 <= 0 || zq <= z0 || zu <= z0 || za <= z0) {
//		usrerr ("height not positive; z0=%f\tzq=%f\tzu=%\tza=%f",
//		       z0, zq, zu, za);
		ier = -2;
		return (ier);
	}

	/* temperatures are Kelvin */
	if (ta <= 0 || ts <= 0) {
//		usrerr ("temps not K; ta=%f\tts=%f", ta, ts);
		ier = -2;
		return (ier);
	}

	/* pressures must be positive */
	if (ea <= 0 || es <= 0 || press <= 0 || ea >= press || es >= press) {
//		usrerr ("press < 0; ea=%f\tes=%f\tpress=%f", ea, es, press);
		ier = -2;
		return (ier);
	}

	/* vapor pressures can't exceed saturation */
	/* if way off stop */
	if ((es - 25.0) > sati(ts) || (ea - 25.0) > satw(ta)) {
//		usrerr ("vp > sat; es=%f\tessat=%f\tea=%f\teasat=%f",
//			es, sati(ts), ea, sati(ta));
		ier = -2;
		return (ier);
	}
	/* else fix them up */
	if (es > sati(ts)) {
		es = sati(ts);
	}
	if (ea > satw(ta)) {
		ea = satw(ta);
	}

	/*
	 * displacement plane height, eq. 5.3 & 5.4
	 */

	d0 = 2 * PAESCHKE * z0 / 3;

	/*
	 * constant log expressions
	 */

	ltsm = log((zu - d0) / z0);
	ltsh = log((za - d0) / z0);
	ltsv = log((zq - d0) / z0);

	/*
	 * convert vapor pressures to specific humidities
	 */
	qa = SPEC_HUM(ea, press);
	qs = SPEC_HUM(es, press);

	/*
	 * convert temperature to potential temperature
	 */

	ta += DALR * za;

	/*
	 * air density at press, virtual temp of geometric mean
	 * of air and surface
	 */
	
	dens = GAS_DEN(press, MOL_AIR,
			VIR_TEMP(sqrt(ta*ts), sqrt(ea*es), press));

	/*
	 * starting value, assume neutral stability, so psi-functions
	 * are all zero
	 */

	ustar = k * u / ltsm;
	factor = k * ustar * dens;
	*e = (qa - qs) * factor * av / ltsv;
	*h = (ta - ts) * factor * cp * ah / ltsh;

	/*
	 * if not neutral stability, iterate on Obukhov stability
	 * length to find solution
	 */

	iter = 0;
	if (ta != ts) {

		lo = HUGE_VAL;

		do {
			last = lo;

			/*
			 * Eq 4.25, but no minus sign as we define
			 * positive H as toward surface
			 */

			/*
			 * There was an error in the old version of this
			 * line that omitted the cubic power of ustar.
			 * Now, this error has been fixed.
			 */

			lo = ustar * ustar * ustar * dens 
				/ (k * g * (*h/(ta*cp) + 0.61 * *e));

			/*
			 * friction velocity, eq. 4.34'
			 */

			ustar = k * u / (ltsm - psi(zu/lo, SM));

			/*
			 * evaporative flux, eq. 4.33'
			 */

			factor = k * ustar * dens;
			*e = (qa - qs) * factor * av /
				(ltsv - psi(zq/lo, SV));

			/*
			 * sensible heat flus, eq. 4.35'
			 * with sign reversed
			 */

			*h = (ta - ts) * factor * ah * cp /
				(ltsh - psi(za/lo, SH));

			diff = last - lo;

		} while (fabs(diff) > THRESH &&
			fabs(diff/lo) > THRESH &&
			++iter < ITMAX);
	}

	ier = (iter >= ITMAX)? -1 : 0;
	
	xlh = LH_VAP(ts);
	if (ts <= FREEZE)
		xlh += LH_FUS(ts);

	/*
	 * latent heat flux (- away from surf)
	 */
	*le = xlh * *e;

	return (ier);
}
