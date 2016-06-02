//---------------------------------------------------------------------------




#include "_net_rad.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _net_rad -- calculates net allwave radiation
**
** SYNOPSIS
**
**	#include "_snobal.h"
**
**      void
**	_net_rad(void)
**
** DESCRIPTION
**      Calculates net allwave radiation from the net solar radiation
**	incoming thermal/longwave radiation, and the snow surface
**	temperature.
**
** GLOBAL VARIABLES READ
**	I_lw
**	S_n
**	T_s_0
**
** GLOBAL VARIABLES MODIFIED
**	R_n
*/

#include <math.h>

#include "ipw.h"
#include "_snobal.h"
#include "snow.h"
#include "radiation.h"

void
_net_rad(void)
{
	R_n = S_n + (SNOW_EMISSIVITY * (I_lw - STEF_BOLTZ * pow(T_s_0, 4)));
}
