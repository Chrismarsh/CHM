#include	"ipw.h"
#include	"envphys.h"

double
efcon(
	double	k,	/* layer thermal conductivity (J/(m K sec)) */
	double	t,	/* layer temperature (K)		    */
	double	p)	/* air pressure (Pa)  			    */
{
	double	etc;
	double	de;
	double	lh;
	double	e;
	double	q;

	/*	calculate effective layer diffusion
		(see Anderson, 1976, pg. 32)		*/
	de = DIFFUS(p, t);

	/*	set latent heat from layer temp.	*/
	if(t > FREEZE)
		lh = LH_VAP(t);
	else if(t == FREEZE)
		lh = (LH_VAP(t) + LH_SUB(t)) / 2.0;
	else
		lh = LH_SUB(t);

	/*	set mixing ratio from layer temp.	*/
	e = sati(t);
	q = MIX_RATIO(e, p);

	/*	calculate effective layer conductivity	*/
	etc = k + (lh * de * q);

	return (etc);
}
