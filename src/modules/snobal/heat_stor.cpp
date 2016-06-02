#include "ipw.h"
#include "envphys.h"

double
heat_stor(
	double	cp,	/* specific heat of layer (J/kg K) */
	double	spm,	/* layer specific mass (kg/m^2)    */
	double	tdif)	/* temperature change (K)          */
{
	double	stor;

	stor = cp * spm * tdif;

	return (stor);
}
