//---------------------------------------------------------------------------




#include "_time_compact.h"
#include "_new_density.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _time_compact -- compact snowcover by gravity over time
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_time_compact(void)
**
** DESCRIPTION
**	This routine "ages" the snowcover by accounting for the compaction
**	or densification by gravity as time passes.  The snowcover's
**	density is increased using the following "half-saturation" function:
**
**		rho(time) = A / (1 + B/time)
**
**	A = "saturation-level" or asymtope which is the maximum density
**	    due to compaction by gravity (approximately 350 kg/m^2)
**	B = the point for half of the saturation level is reached (10 days)
**
**			^
**			|
**		 A: 350 + = = = = = = = = = = = = = = = = = =
**		  	|			*   *
**			|		   *
**	rho		|	       *
**	(kg/m^2)	|	    *
**			|	  *
**	       A/2: 175 + . . . *
**      		|     * .
**      		|   *   .
**      		|  * 	.
**      		| * 	.
**      		|*	.
**      	      0 +-------+----------------------------->
**      		0	B: 10 days		time
**      
**      
** GLOBAL VARIABLES READ
**	time_step
**
** GLOBAL VARIABLES MODIFIED
**	rho
**
*/

#include        "ipw.h"
#include        "_snobal.h"

#define	A	350
	/*
	 *  Maximum density due to compaction by gravity (kg/m^2).
	 */

#define	B	864000
	/*
	 *  Time when half "saturation", i.e., maximum density is reached
	 *  (seconds).
	 *  (864000 = 10 days * 24 hours/day * 60 mins/hr * 60 secs/min)
	 */


void
_time_compact(void)
{	
	double	time;	/* point on time axis corresponding to current
			   density */


	/*
	 *  If the snow is already at or above the maximum density due
	 *  compaction by gravity, then just leave.
	 */
	if ((!snowcover) || (rho > A))
		return;

	/*
	 *  Given the current density, determine where on the time axis
	 *  we are (i.e., solve the function above for "time").
	 */
	time = B / (A / rho - 1);

	/*
	 *  Move along the time axis by the time step, and calculate the
	 *  density at this new time.
	 */
	rho = A / (1 + B/(time + time_step));

        /*
	 *  Adjust the snowcover for this new density.
	 */
	_new_density();
}
