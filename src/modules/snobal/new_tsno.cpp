#include "ipw.h"
#include "snow.h"

double
new_tsno(
	double	spm,	/* layer's specific mass (kg/m^2) 	 */
	double	t0,	/* layer's last temperature (K) 	 */
	double	ccon)	/* layer's adjusted cold content (J/m^2) */
{
	double	tsno;
	double	cp;
	double	tdif;

	cp = CP_ICE(t0);

	tdif = ccon / (spm * cp);
	tsno = tdif + FREEZE;

	return (tsno);
}
