//---------------------------------------------------------------------------




#include "_adj_snow.h"
#include "_adj_layers.h"
#include "_layer_mass.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _adj_snow -- adjust the snowcover with changes in depth and/or mass
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_adj_snow(
**	    double delta_z_s,	|* change in snowcover's depth *|
**	    double delta_m_s)	|* change is snowcover's mass *|
**
** DESCRIPTION
**      This routine adjusts the snowcover for a change in its depth or
**	its mass or both.  The snowcover's density is updated.  If there
**	is a change in the snowcover's depth, the # of layers is recomputed.
**	If there's just a change in the snowcover's mass with no change in
**	its depth, then just the specific masses for the layers are updated.
**
**	The routine ensures that the snowcover's density does NOT exceed
**	a maximum density (currently 750 kg/m^3).  If the adjustments to
**	the snowcover, for some reason, lead to an excessive density, the
**	density is clipped at the maximum, and the depth re-adjusted
**	accordingly.
**
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**	m_s	
**	rho
**	z_s
**
**	(and those variables modified by "_adj_layers" and "_layer_mass")
*/

#include "ipw.h"
#include "_snobal.h"
#include "envphys.h"

#define MAX_SNOW_DENSITY	750
/*
 *  Maximum snow density (kg/m^3)
 */


void
_adj_snow(
	double	delta_z_s,	/* change in snowcover's depth */
	double	delta_m_s)	/* change is snowcover's mass */
{
	/*
	 *  Update depth, mass, and then recompute density.
	 */
	z_s += delta_z_s;
	m_s += delta_m_s;

	if (z_s != 0.0) {
		rho = m_s / z_s;
	} else {
		rho = 0;
	}

	/*
	 *  Clip density at maxium density if necessary.
	 */
	if (rho > MAX_SNOW_DENSITY)
		{
		rho = MAX_SNOW_DENSITY;
		z_s = m_s / rho;
		_adj_layers();
		}
	else	
		{
		/*
		 *  If a change in depth, adjust the layers' depths and masses.
		 */
		if (delta_z_s != 0.0)
			_adj_layers();
		else
		/*
	 	 *  Just a change in the snowcover's mass, so update the
	 	 *  layers' masses.
	 	 */
			_layer_mass();
		}
}
