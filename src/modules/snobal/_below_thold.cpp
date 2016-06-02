//---------------------------------------------------------------------------




#include "_below_thold.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _below_thold -- is a layer's mass below a threshold ?
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      int
**	_below_thold(
**	    double  threshold)	|* current timestep's threshold for a 
**				   layer's mass *|
**
** DESCRIPTION
**      This routine determines if any individual layer's mass is below
**	a given threshold for the current timestep.
**
** RETURN VALUE
**      	1	A layer's mass is less than the threshold.
**
**		0	All layers' masses are greater than the threshold.
**
** GLOBAL VARIABLES READ
**	layer_count
**	m_s
**	m_s_0
**	m_s_l
**
** GLOBAL VARIABLES MODIFIED
*/

#include        "ipw.h"
#include        "_snobal.h"
#include        "envphys.h"

int
_below_thold(
	double	threshold)	/* current timestep's threshold for a
				   layer's mass */
{
	if (layer_count == 0)
		return 0;
	if (layer_count == 1)
		return (m_s < threshold);
	else  /* layer_count == 2 */
		return ((m_s_0 < threshold) || (m_s_l < threshold));
}
