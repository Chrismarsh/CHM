//01/26/07
//---------------------------------------------------------------------------



#include "init_snow.h"
#include "_calc_layers.h"
#include "_layer_mass.h"
#include "_cold_content.h"
#include "_snobal.h"
#include "snow.h"

//---------------------------------------------------------------------------



/*
** NAME
**      init_snow -- initialize the properties for the snowcover
**
** SYNOPSIS
**      #include "snobal.h"
**
**      void
**	init_snow(void)
**
** DESCRIPTION
**      This routine initializes the properties for the snowcover.  It
**	determines the number of layers, their individual properties,
**	the cold content for the snowcover and its layers, and the
**	snowcover's water content.  The following global variables
**	should be initialized before invoking this routine:
**
**		z_s	depth of snowcover (m)
**		rho	density of snowcover (kg/m^3)
**		T_s	average temperature of snowcover (K)
**		T_s_0	temperature of surface layer of snowcover (K)
**		T_s_l	temperature of lower layer of snowcover (K)
**		h2o_sat	% of liquid h2o saturation (0 to 1.0)
**
**		max_h2o_vol	maximum liquid h2o content as volume ratio:
**				    V_water/(V_snow - V_ice) (unitless)
**
** GLOBAL VARIABLES READ
**	h2o_sat
**	layer_count
**	m_s_0
**	m_s_l
**	max_h2o_vol
**	rho
**	T_s
**	T_s_0
**	T_s_l
**	z_s
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
*/

void
init_snow(void)
{
	double	rho_dry;	/* snow density without H2O */

	m_s = rho * z_s;

	_calc_layers();

	if (layer_count == 0) {
		/*
		 *  If mass > 0, then it must be below threshold.
		 *  So turn this little bit of mass into water.
		 */
		if (m_s > 0.0)
			h2o_total += m_s;

		rho = 0.0;
		m_s   = cc_s   = 0.0;
		m_s_0 = cc_s_0 = 0.0;
		m_s_l = cc_s_l = 0.0;

		/*
		 *  Note: Snow temperatures are set to MIN_SNOW_TEMP
		 *	  (as degrees K) instead of 0 K to keep quantization
		 *	  range in output image smaller.
		 */
		T_s = T_s_0 = T_s_l = MIN_SNOW_TEMP + FREEZE;
		h2o_vol = h2o = h2o_max = h2o_sat = 0.0;
	}

	else {
		/*
		 *  Compute specific mass for each layer.
		 */
		_layer_mass();

		cc_s_0 = _cold_content(T_s_0, m_s_0);

		if (layer_count == 2) {
			cc_s_l = _cold_content(T_s_l, m_s_l);
		}
		else {
			T_s_l = MIN_SNOW_TEMP + FREEZE;
			cc_s_l = 0.0;
		}

		/*
		 *  Compute liquid water content as volume ratio, and
		 *  snow density without water.
		 */
		h2o_vol = h2o_sat * max_h2o_vol;
		rho_dry = DRY_SNO_RHO(rho, h2o_vol);

		/*
		 *  Determine maximum liquid water content (as specific mass)
		 *  and the actual liquid water content (as specific mass).
		 */
		h2o_max = H2O_LEFT(z_s, rho_dry, max_h2o_vol);
		h2o = h2o_sat * h2o_max;
	}
}
