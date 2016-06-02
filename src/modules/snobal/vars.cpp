//01/26/07
//---------------------------------------------------------------------------



#include "vars.h"
#include "_snobal.h"

//---------------------------------------------------------------------------



/*
** NAME
**      vars.c - public global variables for snow-balance library
**
** DESCRIPTION
**      Defines and allocates space for variables declared in "snobal.h"
**
*/

#include        "snobal.h"

/* --------------------------------------------------------------------- */

/*
 *  Global variables that are used to communicate with the snobal library
 *  routines.  These are public, i.e., these can be accessed from outside
 *  the library.
 */

/*   variables that control model execution   */

	int	run_no_snow;	/* continue model even if snow disappears? */
	int	stop_no_snow;	/* stopped model because no snow left? */

	void	(*out_func)(void);	/* -> output function */


/*   constant model parameters  */

	double  max_z_s_0;      /* maximum active layer thickness (m) */
	double  max_h2o_vol;    /* max liquid h2o content as volume ratio:
				     V_water/(V_snow - V_ice) (unitless) */


/*   time step information */

	TSTEP_REC  tstep_info[4]; 	/* array of info for each timestep:
						   0 : data timestep
						   1 : normal run timestep
						   2 : medium  "     "
						   3 : small   "     "
					 */

	double	time_step;	/* length current timestep (sec) */
	double  current_time;   /* start time of current time step (sec) */
	double	time_since_out;	/* time since last output record (sec) */


/*   snowpack information   */

	int     layer_count;    /* number of layers in snowcover: 0, 1, or 2 */
	double  z_s;            /* total snowcover thickness (m) */
	double  z_s_0;          /* active layer depth (m) */
	double  z_s_l;          /* lower layer depth (m) */
	double  rho;            /* average snowcover density (kg/m^3) */
	double  m_s;            /* snowcover's specific mass (kg/m^2) */
	double  m_s_0;          /* active layer specific mass (kg/m^2) */
	double  m_s_l;          /* lower layer specific mass (kg/m^2) */
	double  T_s;            /* average snowcover temp (K) */
	double  T_s_0;          /* active snow layer temp (K) */
	double  T_s_l;          /* lower layer temp (C) */
	double  cc_s;           /* snowcover's cold content (J/m^2) */
	double  cc_s_0;         /* active layer cold content (J/m^2) */
	double  cc_s_l;         /* lower layer cold content (J/m^2) */
	double  h2o_sat;        /* % of liquid H2O saturation (relative water
				     content, i.e., ratio of water in snowcover
				     to water that snowcover could hold at
				     saturation) */
	double  h2o_vol;        /* liquid h2o content as volume ratio:
				     V_water/(V_snow - V_ice) (unitless) */
	double  h2o;            /* liquid h2o content as specific mass
				     (kg/m^2) */
	double  h2o_max;        /* max liquid h2o content as specific mass
				     (kg/m^2) */
	double  h2o_total;      /* total liquid h2o: includes h2o in snowcover,
				     melt, and rainfall (kg/m^2) */


/*   climate-data input records   */

	int     ro_data;        /* runoff data? */

	INPUT_REC  input_rec1;	/* input data for start of data timestep */
	INPUT_REC  input_rec2;	/*   "     "   "  end   "   "      "     */

/*   climate-data input values for the current run timestep */

	double  S_n;		/* net solar radiation (W/m^2) */
	double  I_lw;           /* incoming longwave (thermal) rad (W/m^2) */
	double  T_a;            /* air temp (C) */
	double  e_a;            /* vapor pressure (Pa) */
	double  u;              /* wind speed (m/sec) */
	double  T_g;            /* soil temp at depth z_g (C) */
	double  ro;             /* measured runoff (m/sec) */


/*   other climate input   */

	double  P_a;            /* air pressure (Pa) */


/*   measurement heights/depths   */

	int	relative_hts;	/* TRUE if measurements heights, z_T
				   and z_u, are relative to snow
				   surface; FALSE if they are
				   absolute heights above the ground */
	double  z_g;            /* depth of soil temp meas (m) */
	double  z_u;            /* height of wind measurement (m) */
	double  z_T;            /* height of air temp & vapor pressure
				   measurement (m) */
	double  z_0;            /* roughness length */


/*   precipitation info for the current DATA timestep    */

	int	precip_now;	/* precipitation occur for current timestep? */
	double  m_pp;		/* specific mass of total precip (kg/m^2) */
	double  percent_snow;	/* % of total mass that's snow (0 to 1.0) */
	double  rho_snow;       /* density of snowfall (kg/m^3) */
	double  T_pp;           /* precip temp (C) */
	double	T_rain;		/* rain's temp (K) */
	double	T_snow;		/* snowfall's temp (K) */
	double  h2o_sat_snow;   /* snowfall's % of liquid H2O saturation */

/*   precipitation info adjusted for current run timestep   */

	double	m_precip;	/* specific mass of total precip (kg/m^2) */
	double	m_rain;		/*    "      "   of rain in precip (kg/m^2) */
	double	m_snow;		/*    "      "   "  snow "    "    (kg/m^2) */
	double	z_snow;		/* depth of snow in precip (m) */


/*   energy balance info for current timestep        */

	double  R_n;            /* net allwave radiation (W/m^2) */
	double  H;              /* sensible heat xfr (W/m^2) */
	double  L_v_E;          /* latent heat xfr (W/m^2) */
	double  G;              /* heat xfr by conduction & diffusion from soil
				     to snowcover (W/m^2) */
	double  G_0;            /* heat xfr by conduction & diffusion from soil
				     or lower layer to active layer (W/m^2) */
	double  M;              /* advected heat from precip (W/m^2) */
	double  delta_Q;        /* change in snowcover's energy (W/m^2) */
	double  delta_Q_0;      /* change in active layer's energy (W/m^2) */

/*   averages of energy balance vars since last output record   */

	double	R_n_bar;
	double	H_bar;
	double	L_v_E_bar;
	double	G_bar;
	double	G_0_bar;
	double	M_bar;
	double	delta_Q_bar;
	double	delta_Q_0_bar;


/*   mass balance vars for current timestep        */

	double  melt;       	/* specific melt (kg/m^2 or m) */
	double  E;		/* mass flux by evap into air from active
				     layer (kg/m^2/s) */
	double  E_s;		/* mass of evap into air & soil from snowcover
				     (kg/m^2) */
	double  ro_predict;     /* predicted specific runoff (m/sec) */

/*   sums of mass balance vars since last output record   */

	double	melt_sum;
	double	E_s_sum;
	double	ro_pred_sum;
 