//---------------------------------------------------------------------------




#include "_new_density.h"
#include "_adj_layers.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _new_density -- adjust the snowcover's depth and layers for new density
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_new_density(void)
**
** DESCRIPTION
**      This routine adjusts the snowcover's depth for a new density.  The
**	layers are also adjusted accordingly.
**
** GLOBAL VARIABLES READ
**	m_s	
**	rho
**
** GLOBAL VARIABLES MODIFIED
**	z_s
**
**	(and those variables modified by "_adj_layers")
*/

#include "ipw.h"
#include "_snobal.h"

void
_new_density(void)
{
	z_s = m_s / rho;

	_adj_layers();
}
