//---------------------------------------------------------------------------




#include "_h2o_compact.h"
#include "_new_density.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _h2o_compact -- compact snowcover due to liquid H2O that was added
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_h2o_compact(void)
**
** DESCRIPTION
**	This routine compacts or densifies the snowcover based on the
**	amount of liquid H2O that was added to the snowcover from melting
**	and rain.  The snowcover's density is increased using the
**	following "half-saturation" function:
**
**		delta_rho(h2o_added) = A / (1 + B/h2o_added)
**
**	A = "saturation-level" or asymtope which is the difference between
**	    the maximum density due to compaction by liquid H2O
**	    (approximately 550 kg/m^2) and the current density
**	B = the point for half of the saturation level is reached (5 %)
**	    (h2o_added = ratio of mass of liquid h2o added by melting and
**		         rain to the mass of the snowcover)
**
**			^
**			|
**		      A + = = = = = = = = = = = = = = = = = =
**	(550 - current  |			*   *
**	       density)	|		   *
**			|	       *
**	delta_rho	|	    *
**	(kg/m^2)	|	  *
**	            A/2 + . . . *
**      		|     * .
**      		|   *   .
**      		|  * 	.
**      		| * 	.
**      		|*	.
**      	      0 +-------+-----------------------------+	  h2o_added
**      		0	B: 5 %			     1.0
**      
**      
** GLOBAL VARIABLES READ
**	m_rain
**	m_s
**	melt
**
** GLOBAL VARIABLES MODIFIED
**	rho
**
*/

#include        "ipw.h"
#include        "_snobal.h"

#define	MAX_DENSITY	550
	/*
	 *  Maximum density due to compaction by liquid H2O added (kg/m^2).
	 */

#define	B	0.05
	/*
	 *  ratio where half the difference between maximum density and
	 *  current density is reached (ratio from 0.0 to 1.0).
	 */


void
_h2o_compact(void)
{	
	double  A;		/* difference between maximum & current
				   densities */
	double	h2o_added;	/* ratio of mass of liquid H2O added from
			   	   melting and rain to mass of snowcover */

	/*
	 *  If the snow is already at or above the maximum density due
	 *  compaction by liquid H2O, then just leave.
	 */
	if ((!snowcover) || (rho > MAX_DENSITY))
		return;

	A = MAX_DENSITY - rho;
	if (precip_now)
		h2o_added = (melt + m_rain) / m_s;
	else
		h2o_added = melt / m_s;
	if (h2o_added > 0.000001) {
		rho += A / (1 + B/h2o_added);

	        /*
		 *  Adjust the snowcover for this new density.
		 */
		_new_density();
	}
}
