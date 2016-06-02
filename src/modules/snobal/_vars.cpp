//01/26/07
//---------------------------------------------------------------------------



#include "_snobal.h"
#include "_vars.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _vars.c - private global variables for snow-balance library
**
** DESCRIPTION
**      Defines and allocates space for variables declared in "_snobal.h"
**
*/

/* --------------------------------------------------------------------- */

/*
 *  Global variables that are private (internal) to the snobal library.
 */

	INPUT_REC input_deltas[4];	/* deltas for climate-input parameters
					   over each timestep */

	PRECIP_REC precip_info[4];	/* array of precip info adjusted for
					   each timestep */

	int  computed[4];		/* array of flags for each timestep;
					   TRUE if computed values for input
					   deltas and precip arrays */

	int	isothermal;	/* melting? */
        int     snowcover;      /* snow on gnd at start of current timestep? */
 