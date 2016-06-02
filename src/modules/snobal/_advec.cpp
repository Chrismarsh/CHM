//---------------------------------------------------------------------------




#include "_advec.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _advec -- calculates advected energy at a point
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_advec(void)
**
** DESCRIPTION
**      This routine calculates the advected energy for a 2-layer snowcover
**	if there's precipitation for the current timestep.
**      
** GLOBAL VARIABLES READ
**	m_rain
**	m_snow
**	precip_now
**	T_rain
**	T_s_0
**	T_snow
**	time_step
**
** GLOBAL VARIABLES MODIFIED
**	M
*/

#include        "ipw.h"
#include        "_snobal.h"
#include	"envphys.h"

void
_advec(void)
{
	if (precip_now) {
		M = (heat_stor(CP_WATER(T_rain), m_rain, (T_rain - T_s_0)) +
		     heat_stor(CP_ICE(T_snow), m_snow, (T_snow - T_s_0)))
 		    / time_step;
	}
	else
		M = 0.0;
}
