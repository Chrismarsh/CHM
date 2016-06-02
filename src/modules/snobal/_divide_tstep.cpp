//---------------------------------------------------------------------------




#include "_divide_tstep.h"
#include "_snobal.h"
#include "_below_thold.h"
#include "_do_tstep.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _divide_tstep -- divide a timestep into smaller timesteps
**
** SYNOPSIS
**	#include "_snobal.h"
**
**	int
**	_divide_tstep(
**	    TSTEP_REC *tstep;	|* record of timestep to be divided *|
**
** DESCRIPTION
**	This routine divides a timestep into smaller timesteps.  For
**	each of these smaller timestep, the routine either run the
**	model for that timestep, or further subdivides that timestep
**	into even smaller timesteps.
**
**	The routine will set the flag 'stop_no_snow' to TRUE if
**
**		a)  the output function pointed to by 'out_func' is called, and
**		b)  the flag 'run_no_snow' is FALSE, and
**		c)  there is no snow remaining on the ground at the end of
**		    timestep
**
** RETURN VALUE
**
**	TRUE	The timestep was successfully divided into smaller timesteps.
**
**	FALSE	An error occured while running the model during one of the
**		smaller timesteps.  An message explaining the error has
**		been stored with the 'usrerr' routine.
**
** GLOBAL VARIABLES READ 
**	layer_count
**	precip_now
**	ro_data
**	tstep_info
**
** GLOBAL VARIABLES MODIFIED
*/

#include	"ipw.h"
#include        "_snobal.h"

int
_divide_tstep(
	TSTEP_REC *tstep)	/* record of timestep to be divided */
{
	int	    next_level; 	/* # of next level of timestep */
	TSTEP_REC  *next_lvl_tstep;	/* info of next level of timestep */
	INPUT_REC  *curr_lvl_deltas;	/* -> input-deltas of current level */
	INPUT_REC  *next_lvl_deltas;	/* -> input-deltas of next level */
	PRECIP_REC *curr_lvl_precip;	/* -> precip data of current level */
	PRECIP_REC *next_lvl_precip;	/* -> precip data of next level */
	int	    i;			/* loop index */


	/*
	 *  Fetch the record for the timestep at the next level.
	 */
	next_level = tstep->level + 1;
	next_lvl_tstep = tstep_info+ next_level;

	curr_lvl_deltas = input_deltas + tstep->level;
	next_lvl_deltas = input_deltas + next_level;

	curr_lvl_precip = precip_info + tstep->level;
	next_lvl_precip = precip_info + next_level;

	/*
	 *  If this is the first time this new level has been used during
	 *  the current data timestep, then calculate its input deltas
	 *  and precipitation values.
	 */
	if (! computed[next_level]) {
		next_lvl_deltas->S_n  = curr_lvl_deltas->S_n /
						    next_lvl_tstep->intervals;
		next_lvl_deltas->I_lw = curr_lvl_deltas->I_lw /
						    next_lvl_tstep->intervals;
		next_lvl_deltas->T_a  = curr_lvl_deltas->T_a /
						    next_lvl_tstep->intervals;
		next_lvl_deltas->e_a  = curr_lvl_deltas->e_a /
						    next_lvl_tstep->intervals;
		next_lvl_deltas->u    = curr_lvl_deltas->u /
						    next_lvl_tstep->intervals;
		next_lvl_deltas->T_g  = curr_lvl_deltas->T_g /
						    next_lvl_tstep->intervals;
		if (ro_data)
			next_lvl_deltas->ro = curr_lvl_deltas->ro /
						    next_lvl_tstep->intervals;

		if (precip_now) {
			next_lvl_precip->m_pp   = curr_lvl_precip->m_pp /
						    next_lvl_tstep->intervals;
			next_lvl_precip->m_rain = curr_lvl_precip->m_rain /
						    next_lvl_tstep->intervals;
			next_lvl_precip->m_snow = curr_lvl_precip->m_snow /
						    next_lvl_tstep->intervals;
			next_lvl_precip->z_snow = curr_lvl_precip->z_snow /
						    next_lvl_tstep->intervals;
		}

		computed[next_level] = 1;
	}

	/*
	 *  For each the new smaller timestep, either subdivide them if
	 *  below their mass threshold, or run the model for them.
	 */
	for (i = 0; (i < next_lvl_tstep->intervals) && !stop_no_snow; i++)
	    if ((next_level != SMALL_TSTEP) &&
				_below_thold(next_lvl_tstep->threshold)) {
		if (! _divide_tstep(next_lvl_tstep))
			return 0;
	    }
	    else {
		if (! _do_tstep(next_lvl_tstep))
			return 0;
	    }

	/*
	 *  Output if this timestep is divided?
	 */
	if (tstep->output & DIVIDED_TSTEP) {
	  	(*out_func)();
		if (!run_no_snow && (layer_count == 0))
			stop_no_snow = 1;
	}

	return 1;
}
