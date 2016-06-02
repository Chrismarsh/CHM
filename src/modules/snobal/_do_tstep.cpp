//---------------------------------------------------------------------------




#include "_do_tstep.h"
#include "_snobal.h"
#include "_e_bal.h"
#include "_mass_bal.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _do_tstep -- do model calculations for 1 timestep
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      int
**	_do_tstep(
**	    TSTEP_REC *tstep;  |* timestep's record *|
**
** DESCRIPTION
**      This routine performs the model's calculations for a single timestep.
**      It requires that these climate variables have been initialized:
**      
**      	S_n
**      	I_lw
**      	T_a
**      	e_a
**      	u
**      	T_g
**      	ro
**
**	The routine also requires the precipitation data have been adjusted
**	for the timestep, and have been stored in the array:
**
**		precip_info
**
**	if the flag 'precip_now' is TRUE.  The routine will set the flag
**	'stop_no_snow' to TRUE if
**
**		a)  the output function pointed to by 'out_func' is called, and
**		b)  the flag 'run_no_snow' is FALSE, and
**		c)  there is no snow remaining on the ground at the end of
**		    timestep
**
** RETURN VALUE
**
**	TRUE	The model's calculations were completed.
**
**	FALSE	An error occured, and a message explaining the error has
**		been stored with the 'usrerr' routine.
**
** GLOBAL VARIABLES READ
**	delta_Q
**	delta_Q_0
**	E_s
**	G
**	H
**	L_v_E
**	layer_count
**	M
**	melt
**	out_func
**	precip_now
**	R_n
**	ro_data
**	ro_predict
**	run_no_snow
**
** GLOBAL VARIABLES MODIFIED
**	current_time
**	curr_time_hrs
**	delta_Q_0_bar
**	delta_Q_bar
**	e_a
**	E_s_sum
**	G_bar
**	H_bar
**	h2o_total
**	I_lw
**	L_v_E_bar
**	M_bar
**	m_precip
**	m_rain
**	m_snow
**	melt_sum
**	R_n_bar
**	ro
**	ro_pred_sum
**	S_n
**	snowcover
**	T_a
**	T_g
**	u
**	time_since_out
**	time_step
**	z_snow
*/

#include        "ipw.h"
#include        "_snobal.h"

/*
 *  A macro to update a time-weighted average for a quantity.
 *	avg		current average
 *	total_time	the time interval the current average applies to
 *	value		new value to be averaged in
 *	time_incr	the time interval the new value applies to
 */
#define TIME_AVG(avg,total_time,value,time_incr) \
		( ((avg) * (total_time) + (value) * (time_incr)) \
		/ ((total_time) + (time_incr)) )

int
_do_tstep(
	TSTEP_REC *tstep)  /* timestep's record */
{
	time_step = tstep->time_step;

	if (precip_now) {
		m_precip = precip_info[tstep->level].m_pp;
		m_rain   = precip_info[tstep->level].m_rain;
		m_snow   = precip_info[tstep->level].m_snow;
		z_snow   = precip_info[tstep->level].z_snow;
	}

	h2o_total = 0.0;

	/*
	 *  Is there a snowcover?
	 */
	snowcover = (layer_count > 0);

	/*
	 *  Calculate energy transfer terms
	 */
	if (! _e_bal())
		return 0;

	/*
	 *  Adjust mass and calculate runoff
	 */
	_mass_bal();

	/*
	 *  Update the averages for the energy terms and the totals for mass
	 *  changes since the last output.
	 */
	if (time_since_out > 0.0) {
		R_n_bar       = TIME_AVG(R_n_bar, 	time_since_out,
					 R_n, 		time_step);
		H_bar         = TIME_AVG(H_bar, 	time_since_out,
					 H, 		time_step);
		L_v_E_bar     = TIME_AVG(L_v_E_bar, 	time_since_out,
					 L_v_E, 	time_step);
		G_bar         = TIME_AVG(G_bar, 	time_since_out,
					 G, 		time_step);
		M_bar         = TIME_AVG(M_bar, 	time_since_out,
					 M, 		time_step);
		delta_Q_bar   = TIME_AVG(delta_Q_bar,	time_since_out,
					 delta_Q, 	time_step);
		G_0_bar       = TIME_AVG(G_0_bar, 	time_since_out,
					 G_0, 		time_step);
		delta_Q_0_bar = TIME_AVG(delta_Q_0_bar, time_since_out,
					 delta_Q_0, 	time_step);

		E_s_sum     += E_s;
		melt_sum    += melt;
		ro_pred_sum += ro_predict;

		time_since_out += time_step;
	}
	else {
		R_n_bar = R_n;
		H_bar = H;
		L_v_E_bar = L_v_E;
		G_bar = G;
		M_bar = M;
		delta_Q_bar = delta_Q;
		G_0_bar = G_0;
		delta_Q_0_bar = delta_Q_0;

		E_s_sum     = E_s;
		melt_sum    = melt;
		ro_pred_sum = ro_predict;

		time_since_out = time_step;
	}

	/* increment time */
	current_time += time_step;

	if (tstep->output & WHOLE_TSTEP) {
		(*out_func)();
		if (!run_no_snow && (layer_count == 0))
			stop_no_snow = 1;
	}

	/*
	 *  Update the model's input parameters
	 */
	S_n  += input_deltas[tstep->level].S_n;
	I_lw += input_deltas[tstep->level].I_lw;
	T_a  += input_deltas[tstep->level].T_a;
	e_a  += input_deltas[tstep->level].e_a;
	u    += input_deltas[tstep->level].u;
	T_g  += input_deltas[tstep->level].T_g;
	if (ro_data)
		ro += input_deltas[tstep->level].ro;

	return 1;
}
