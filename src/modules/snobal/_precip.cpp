//---------------------------------------------------------------------------




#include "_precip.h"
#include "init_snow.h"
#include "_adj_snow.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _precip -- process a precipitation event
** 
** SYNOPSIS
**	#include "_snobal.h"
**
**      void
**	_precip(void)
** 
** DESCRIPTION
**      This routine processes a precipitation event, i.e., the current
**	precip record, if there's one for the current timestep.  It
**	determines if the precip is rain or snow which increases the
**	snowcover.
** 
** GLOBAL VARIABLES READ
**	h2o_sat_snow
**	m_rain
**	m_precip
**	max_h2o_vol
**	precip_now
**	rho_snow
**	snowcover
**	T_snow
**	z_snow
**
** GLOBAL VARIABLES MODIFIED
**	h2o
**	h2o_sat
**	h2o_total
**	rho
**	T_s
**	T_s_0
**	z_s
*/

#include "ipw.h"
#include "snow.h"
#include "_snobal.h"

void
_precip(void)
{
	double	h2o_vol_snow;	/* liquid water content of new snowfall as
				   volume ratio */

	if (precip_now) {
		if (snowcover) {
			/*
			 *  Adjust snowcover's depth and mass by snowfall's
			 *  depth and the total precipitation mass.
			 */
			_adj_snow(z_snow, m_precip);

			/*
			 *  Determine the additional liquid water that's in
			 *  the snowfall, and then add its mass to liquid
			 *  water in the whole snowcover.
			 */
			h2o_vol_snow = h2o_sat_snow * max_h2o_vol;
			h2o += H2O_LEFT(z_snow, rho_snow, h2o_vol_snow);
		}
                else {
			/*
			 *  Use snowfall, if any, to setup a new snowcover.
			 */
			if (m_snow > 0.0) {
				z_s = z_snow;
				rho = rho_snow;
				T_s = T_snow;
				T_s_0 = T_snow;
				T_s_l = T_snow;
				h2o_sat = h2o_sat_snow;

				init_snow();
			}
		}

		/*
		 *  Add rainfall and water in the snowcover to total
		 *  liquid water.
		 */
		h2o_total += h2o + m_rain;
        }
        else
		/*
		 *  Add water in the snowcover to total liquid water.
		 */
		h2o_total += h2o;
}
