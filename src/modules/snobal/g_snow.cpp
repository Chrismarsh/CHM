//---------------------------------------------------------------------------




#include "g_snow.h"

//---------------------------------------------------------------------------



#include "ipw.h"
#include "snow.h"

double
g_snow(
	double	rho1,	/* upper snow layer's density (kg/m^3)	*/
	double	rho2,	/* lower  "     "        "    (kg/m^3)	*/
	double	ts1,	/* upper snow layer's temperature (K)	*/
	double	ts2,	/* lower  "     "         "       (K)	*/
	double	ds1,	/* upper snow layer's thickness (m)	*/
	double	ds2,	/* lower  "     "         "     (m)	*/
	double	pa)	/* air pressure (Pa)			*/
{
	double	kcs1;
	double	kcs2;
	double	k_s1;
	double	k_s2;
	double	g;


/*	calculate G	*/
	if (ts1 == ts2)
		g = 0.0;
	else {
	/*	set snow conductivity	*/
		kcs1 = KTS(rho1);
		kcs2 = KTS(rho2);
		k_s1 = efcon(kcs1, ts1, pa);
		k_s2 = efcon(kcs2, ts2, pa);

	/*	calculate g	*/
		g = ssxfr(k_s1, k_s2, ts1, ts2, ds1, ds2);
	}

	return (g);
}
