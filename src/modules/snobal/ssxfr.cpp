#include "ipw.h"
#include "envphys.h"

double
ssxfr(
	double	k1,	/* layer 1's thermal conductivity (J / (m K sec))  */
	double	k2,	/* layer 2's    "         "                        */
	double	t1,	/* layer 1's average layer temperature (K)	   */
	double	t2,	/* layer 2's    "      "        "         	   */
	double	d1,     /* layer 1's thickness (m)			   */
	double	d2)     /* layer 2's    "       "			   */
{
	double	xfr;

	xfr = 2.0 * (k1 * k2 * (t2 - t1)) / ((k2 * d1) + (k1 * d2));

	return (xfr);
}
