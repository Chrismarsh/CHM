//---------------------------------------------------------------------------




#include "_cold_content.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _cold_content -- calculates cold content for a layer
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      double
**	_cold_content(
**	    double  temp,		|* temperature of layer *|
**	    double  mass)		|* specific mass of layer *|
**
** DESCRIPTION
**      This routine calculates the cold content for a layer (i.e., the
**	energy required to bring its temperature to freezing) from the
**	layer's temperature and specific mass.
** 
** RETURN VALUE
**	The layer's cold content.
*/

#include        "ipw.h"
#include        "_snobal.h"
#include        "envphys.h"

double
_cold_content(
	double	temp,		/* temperature of layer */
	double	mass)		/* specific mass of layer */
{
	if (temp < FREEZE)
	 	return heat_stor(CP_ICE(temp), mass, (temp - FREEZE));
	else
		return 0.0;
}
