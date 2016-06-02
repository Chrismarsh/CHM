//---------------------------------------------------------------------------




#include "_calc_layers.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _calc_layers -- determine # of layers in snowcover and their depths
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_calc_layers(void);
**
** DESCRIPTION
**      This routine determines the # of layers in the snowcover based its
**	depth and mass.  Usually, there are are 2 layers: the surface (active)
**	and the lower layer.  The depth of the surface layer is set to the
**	maximum depth for the surface layer (variable "max_z_s_0").  The
**	remaining depth constitutes the lower layer.  The routine checks
**	to see if the mass of this lower layer is above the minium threshold
**	(i.e., the mass threshold for the small run timestep).  If not,
**	the surface layer is the whole snowcover, and there's no lower
**	layer.
**
** GLOBAL VARIABLES READ
**	m_s
**	max_z_s_0
**	rho
**	tstep_info
**	z_s
**
** GLOBAL VARIABLES MODIFIED
**	layer_count
**	z_s
**	z_s_0
**	z_s_l
*/

#include "ipw.h"
#include "_snobal.h"

void
_calc_layers(void)
{
	if (m_s <= tstep_info[SMALL_TSTEP].threshold) {
		/*
		 *  Less than minimum layer mass, so treat as no snowcover.
		 */
		layer_count = 0;
		z_s = z_s_0 = z_s_l = 0.0;
	}
	else if (z_s < max_z_s_0) {
		/*
		 *  Not enough depth for surface layer and the lower layer,
		 *  so just 1 layer: surface layer.
		 */
		layer_count = 1;
		z_s_0 = z_s;
		z_s_l = 0.0;
	}
	else {
		/*
		 *  Enough depth for both layers.
		 */
		layer_count = 2;
		z_s_0 = max_z_s_0;
		z_s_l = z_s - z_s_0;

		/*
		 *  However, make sure there's enough MASS for the lower
		 *  layer.  If not, then there's only 1 layer.
		 */
		if (z_s_l * rho < tstep_info[SMALL_TSTEP].threshold) {
			layer_count = 1;
			z_s_0 = z_s;
			z_s_l = 0.0;
		}
	}
}
