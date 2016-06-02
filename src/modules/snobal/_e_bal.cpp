//---------------------------------------------------------------------------




#include "_e_bal.h"
#include "_net_rad.h"
#include "_h_le.h"
#include "_advec.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _e_bal -- calculates point energy budget for 2-layer snowcover
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      int
**	_e_bal(void)
**
** DESCRIPTION
**      Calculates point energy budget for 2-layer snowcover.
**      
** RETURN VALUE
**
**	TRUE	The calculations were completed.
**
**	FALSE	An error occured, and a message explaining the error has
**		been stored with the 'usrerr' routine.
**
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**
*/

#include        "ipw.h"
#include        "_snobal.h"
#include        "snow.h"

int
_e_bal(void)
{
	if (snowcover) {

	/**	calculate energy xfer terms  **/

	/*      calculate net radiation */

	        _net_rad();

	/*      calculate H & L_v_E  (and E as well)       */

	        if (! _h_le())
			return 0;

	/*      calculate G & G_0(conduction/diffusion heat xfr)    */

	 	if (layer_count == 1) {
	        	G = g_soil (rho, T_s_0, T_g, z_s_0, z_g, P_a);
	                G_0 = G;
		}
		else {  /*  layer_count == 2  */
	        	G = g_soil (rho, T_s_l, T_g, z_s_l, z_g, P_a);
	                G_0 = g_snow (rho, rho, T_s_0, T_s_l, z_s_0, z_s_l,
					P_a);
		}
	
	/*      calculate advection     */

		_advec();

	/*      sum E.B. terms  */

	        /* surface energy budget */
		delta_Q_0 = R_n + H + L_v_E + G_0 + M;

	        /* total snowpack energy budget */
	        if (layer_count == 1)
	                delta_Q = delta_Q_0;
	        else  /* layer_count == 2 */
	                delta_Q = delta_Q_0 + G - G_0;
	}
	else {
		R_n = 0.0;

		H = L_v_E = E = 0.0;

		G = G_0 = 0.0;

		M = 0.0;

		delta_Q = delta_Q_0 = 0.0;
	}

	return 1;
}
