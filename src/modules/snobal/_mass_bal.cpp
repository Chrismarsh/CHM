//---------------------------------------------------------------------------




#include "_mass_bal.h"
#include "_time_compact.h"
#include "_precip.h"
#include "_snowmelt.h"
#include "_evap_cond.h"
#include "_h2o_compact.h"
#include "_runoff.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _mass_bal -- calculates point mass budget of 2-layer snowcover
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_mass_bal(void)
**
** DESCRIPTION
**      Calculates the point mass budget for 2-layer energy budget snowmelt
**	model.  It then solves for new snow temperatures.
**      
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**
*/

#include        "ipw.h"
#include        "_snobal.h"
#include        "snow.h"

void
_mass_bal(void)
{
        /***    adjust mass and calc. runoff    ***/

	/*	age snow by compacting snow due to time passing */
	_time_compact();

	/*	process precipitation event */
	_precip();

	/*      calculate melt or freezing and adjust cold content */

        _snowmelt();

        /*      calculate evaporation and adjust snowpack       */

	_evap_cond();

	/*	compact snow due to H2O generated (melt & rain) */
	_h2o_compact();

        /*      calculate runoff, and adjust snowcover */

        _runoff();

        /*
	 *  adjust layer temps if there was a snowcover at start of the
	 *  timestep and there's still snow on the ground
	 */
	if (snowcover) {
		if (layer_count == 1) {
	                T_s_0 = new_tsno (m_s_0, T_s_0, cc_s_0);
			T_s = T_s_0;
		}
		else if (layer_count == 2) {
	        	if (isothermal)
	                	T_s = T_s_l = T_s_0 = FREEZE;
	        	else {
	                	T_s_0 = new_tsno (m_s_0, T_s_0, cc_s_0);
	                	T_s_l = new_tsno (m_s_l, T_s_l, cc_s_l);
	                	T_s = new_tsno (m_s, T_s, cc_s);
	        	}
		}
	}
}
