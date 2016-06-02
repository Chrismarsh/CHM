//---------------------------------------------------------------------------




#include "_snowmelt.h"
#include "_adj_snow.h"
#include <cmath>
//---------------------------------------------------------------------------



/*
** NAME
**      _snowmelt -- calculates melt or re-freezing at a point
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_snowmelt(void)
**
** DESCRIPTION
**      Calculates melting or re-freezing for point 2-layer energy balance
**      snowmelt model.
**      
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**
*/

#include        "ipw.h"
#include        "_snobal.h"
#include        "envphys.h"
#include        "snow.h"

void
_snowmelt(void)
{
        double  Q_0;            /* energy available for surface melt */
        double  Q_l;		/* energy available for lower layer melt */
        double  Q_freeze;       /* energy used for re-freezing */
        double  Q_left;         /* energy left after re_freezing */
        double  h2o_refrozen;   /* amount of liquid H2O that was refrozen */


	/*
	 *  If no snow on ground at start of timestep, then just exit.
	 */
	if (!snowcover) {
		melt = 0.0;
		return;
	}

        /***    calculate melt or freezing, and adjust cold content */

        /*** calculate surface melt ***/

	/* energy for surface melt */
        Q_0 = (delta_Q_0 * time_step) + cc_s_0;

        if (Q_0 > 0.0) {
                melt = MELT(Q_0);
                cc_s_0 = 0.0;
        }
        else if (Q_0 == 0.0) {
                melt = 0.0;
                cc_s_0 = 0.0;
        }
        else {
                melt = 0.0;
                cc_s_0 = Q_0;
        }


        /*** calculate lower layer melt ***/

        if (layer_count == 2) {
		/* energy for layer melt */
                Q_l = ((G - G_0) * time_step) + cc_s_l;

	        if (Q_l > 0.0) {
	                melt += MELT(Q_l);
	                cc_s_l= 0.0;
	        }
	        else if (Q_l == 0.0)
	                cc_s_l= 0.0;
	        else
	                cc_s_l= Q_l;
	}
	else {  /* layer_count == 1 */
                Q_l = 0.0;
	}

	h2o_total += melt;


        /*** adjust layers for re-freezing ***/

        /*    adjust surface layer    */

	h2o_refrozen = 0.0;

        if (cc_s_0 < 0.0) {
                /* if liquid h2o present, calc refreezing and adj cc_s_0 */
                if (h2o_total > 0.0) {
                        Q_freeze = h2o_total * (z_s_0/z_s) * LH_FUS(FREEZE);
                        Q_left = Q_0 + Q_freeze;

                        if (Q_left <= 0.0) {
                                h2o_refrozen = h2o_total * (z_s_0/z_s);
                                cc_s_0 = Q_left;
                        }
                        else {
                                h2o_refrozen = (h2o_total * (z_s_0/z_s)) -
								MELT(Q_left);
                                cc_s_0 = 0.0;
                        }
                }
        }

        /*    adjust lower layer for re-freezing */

        if ((layer_count == 2) && (cc_s_l < 0.0)) {
                /* if liquid h2o, calc re-freezing and adj cc_s_l */
                if (h2o_total > 0.0) {
                        Q_freeze = h2o_total * (z_s_l/z_s) * LH_FUS(FREEZE);
                        Q_left = Q_l + Q_freeze;

                        if (Q_left <= 0.0) {
                                h2o_refrozen += h2o_total * (z_s_l/z_s);
                                cc_s_l= Q_left;
                        }
                        else {
                                h2o_refrozen += ((h2o_total* (z_s_l/z_s)) -
								MELT(Q_left));
                                cc_s_l= 0.0;
                        }
                }
        }

	/*
	 *  Note:  because of rounding errors, h2o_refrozen may not
	 * 	   be exactly the same as h2o_total.  Check for this
	 *	   case, and if so, then just zero out h2o_total.
	 */
	if (fabs(h2o_total - h2o_refrozen) <= 1e-8) {
	    h2o_total = 0.0;
	} else {
	    h2o_total -= h2o_refrozen;
	}

	/***	determine if snowcover is isothermal    ***/

	if ((layer_count == 2) && (cc_s_0 == 0.0) && (cc_s_l == 0.0))
		isothermal = 1;
	else if ((layer_count == 1) && (cc_s_0 == 0.0))
		isothermal = 1;
	else
		isothermal = 0;

        /***    adjust depth and density for melt  ***/

        if (melt > 0.0)
		_adj_snow( -(melt/rho), 0.0);

        /***    set total cold content   ***/
        if (layer_count == 2)
                cc_s = cc_s_0 + cc_s_l;
        else if (layer_count == 1)
                cc_s = cc_s_0;
}
