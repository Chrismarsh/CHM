//---------------------------------------------------------------------------




#include "_runoff.h"
#include "_adj_snow.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _runoff -- calculates runoff from snowcover
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_runoff(void)
**
** DESCRIPTION
**      Calculates runoff for point energy budget 2-layer snowmelt model
**      
** GLOBAL VARIABLES READ
**	h2o_total
**	layer_count
**	snowcover
**	max_h2o_vol
**	z_s
**
** GLOBAL VARIABLES MODIFIED
**	h2o
**	h2o_max
**	h2o_sat
**	h2o_vol
**	rho
**	ro_predict
*/

#include        "ipw.h"
#include        "_snobal.h"
#include        "snow.h"

void
_runoff(void)
{
	double	m_s_dry;	/* snowcover's mass without liquid H2O */
	double	rho_dry;	/* snow density without liquid H2O */

        /* calculate runoff */

	/*
	 *  If no snow on ground at start of timestep or no layers currently,
	 *  then all water (e.g., rain) is runoff.
	 */
	if ((!snowcover) || (layer_count == 0)) {
		ro_predict = h2o_total;
		return;
	}

        /*
	 *  Determine the snow density without any water, and the maximum
	 *  liquid water the snow can hold.
	 */
	m_s_dry = m_s - h2o_total;
	rho_dry = m_s_dry / z_s;
	h2o_max = H2O_LEFT(z_s, rho_dry, max_h2o_vol);

	/*
	 *  Determine runoff, and water left in the snow
	 */
        if (h2o_total > h2o_max) {
                ro_predict = h2o_total - h2o_max;
                h2o = h2o_max;
		h2o_sat = 1.0;
		h2o_vol = max_h2o_vol;

		/*
		 *  Update the snowcover's mass for the loss of runoff.
		 */
		_adj_snow(0.0, -ro_predict);
	}
        else {
		ro_predict = 0.0;
		h2o = h2o_total;
		h2o_sat = h2o / h2o_max;
		h2o_vol = h2o_sat * max_h2o_vol;
	}
}
