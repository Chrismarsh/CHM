//---------------------------------------------------------------------------




#include "_layer_mass.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _layer_mass -- calculate the specific mass for each snow layer
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_layer_mass(void)
**
** DESCRIPTION
**      This routine computes the specific mass for each snow layer in
**	the snowcover.  A layer's mass is based its depth and the
**	average snowcover density.
**
** GLOBAL VARIABLES READ
**	layer_count
**	rho
**	z_s_0
**	z_s_l
**
** GLOBAL VARIABLES MODIFIED
**	m_s_0
**	m_s_l
*/

#include "ipw.h"
#include "_snobal.h"

void
_layer_mass(void)
{
	if (layer_count == 0) {
		m_s_0 = 0.0;
		m_s_l = 0.0;
	}
	else {  /* layer_count is 1 or 2 */
		m_s_0 = rho * z_s_0;
		if (layer_count == 2)
			m_s_l = rho * z_s_l;
		else
			m_s_l = 0.0;
	}
}
