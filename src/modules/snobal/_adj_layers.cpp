//---------------------------------------------------------------------------




#include "_adj_layers.h"
#include "_calc_layers.h"
#include "_cold_content.h"
#include "_layer_mass.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _adj_layers -- adjust the layers because of new snowcover depth
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_adj_layers(void);
**
** DESCRIPTION
**      This routine adjusts the layers of the snowcover because the
**	snowcover's depth has changed.  It is assumed that the snowcover's
**	density has already been updated.  The # of layers are recomputed
**	based on the overall snowcover depth.  Their depths and masses
**	are updated as well.  If a layer has been created due to an
**	increase in the snowcover's depth, its temperature and cold content
**	are initialized. 
**
** GLOBAL VARIABLES READ
**	layer_count
**
** GLOBAL VARIABLES MODIFIED
**	cc_s
**	cc_s_0
**	cc_s_l
**	h2o
**	h2o_max
**	h2o_total
**	h2o_vol
**	m_s
**	m_s_0
**	m_s_l
**	rho
**	T_s
**	T_s_0
**	T_s_l
**
**	(and those variables modified by "_calc_layer" and "_layer_mass")
*/

#include "ipw.h"
#include "envphys.h"
#include "_snobal.h"

void
_adj_layers(void)
{
	int prev_layer_count;	/* previous # of layers, if change in depth */

	/*
	 *  Recompute then number of layers and see if there's been
	 *  a change in the # of layers.  Note:  since this routine
	 *  is called to adjust an existing snowcover, the current # of
	 *  layers must be either 1 or 2 while the new # of layers may
	 *  either be 0, 1 or 2.
	 *
	 *	current #	new #
	 *	of layers	of layers
	 *
	 *	   1	   -->	   0
	 *	   1	   -->	   1	(no change)
	 *	   1	   -->	   2
	 *	   2	   -->	   0
	 *	   2	   -->	   1
	 *	   2	   -->	   2	(no change)
	 */
	prev_layer_count = layer_count;  /* must be > 0 */
	_calc_layers();

	if (layer_count == 0) {
		/*
		 *  1 or 2 layers --> 0 layers
		 */
		rho = 0.0;

		/*
		 *  If mass > 0, then it must be below threshold.
		 *  So turn this little bit of mass into water.
		 */
		if (m_s > 0.0)
			h2o_total += m_s;

		m_s   = cc_s   = 0.0;
		m_s_0 = cc_s_0 = 0.0;

		/*
		 *  Note: Snow temperatures are set to MIN_SNOW_TEMP
		 *	  (as degrees K) instead of 0 K to keep quantization
		 *	  range in output image smaller.
		 */
		T_s = T_s_0 = MIN_SNOW_TEMP + FREEZE;

		if (prev_layer_count == 2) {
			m_s_l = cc_s_l = 0.0;
			T_s_l = MIN_SNOW_TEMP + FREEZE;
		}
		h2o_vol = h2o = h2o_max = h2o_sat = 0.0;
	}
 
	else {
		_layer_mass();

		if ((prev_layer_count == 1) && (layer_count == 2)) {
			/*
			 *  1 layer --> 2 layers, add lower layer
			 */
			T_s_l = T_s;
			cc_s_l = _cold_content(T_s_l, m_s_l);
			}

		else if ((prev_layer_count == 2) && (layer_count == 1)) {
			/*
			 *  2 layers --> 1 layer, remove lower layer
			 */
			T_s_l = MIN_SNOW_TEMP + FREEZE;
			cc_s_l = 0.0;
		}
	}

}
