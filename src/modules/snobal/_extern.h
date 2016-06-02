#include "txtio.h"
#include "snobal.h"	/* for typedef "input_rec" */

/***    external var type defs          ***/

/*   file pointers              */

extern  TEXT_FD_T    in;        /* input file ptr */
extern  TEXT_FD_T    snf;       /* snow properties input file ptr */
extern  TEXT_FD_T    mhf;       /* measurement heights input file ptr */
extern  TEXT_FD_T    prf;       /* precipitation input file ptr */
extern  TEXT_FD_T    out;       /* output file ptr */


/*   file names                 */

extern  char    *in_filename;   /* input file name */
extern  char    *sn_filename;   /* snow properties input file name */
extern  char    *mh_filename;   /* measurement heights input file name */
extern  char    *pr_filename;   /* precipitation input file name */
extern  char    *out_filename;  /* output file name */


/*   input constants            */

extern  double  elevation;      /* elevation (m) */
extern  double  max_z_s_0;      /* maximum active layer thickness (m) */
extern  double  max_h2o_vol;    /* max liquid h2o content as volume ratio:
				     V_water/(V_snow - V_ice) (unitless) */
extern	double  min_z;		/* min depth for a snow layer (m) */
extern	double  critical_mass;	/* critical mass threshold for lower layer
				     (kg/m^2) */
extern  int	output_mode;	/* controls how often output records written */


/*   constant model parameters  */

extern  double  P_a;            /* air pressure (Pa) */
extern  double  z_g;            /* depth of soil temp meas (m) */


/*   "in" data file record      */

extern	input_rec  in_rec;


/*   climate data input values	*/

extern  double  S_n;		/* net solar radiation (W/m^2) */
extern  double  I_lw;           /* incoming longwave (thermal) rad (W/m^2) */
extern  double  T_a;            /* air temp (C) */
extern  double  e_a;            /* vapor pressure (Pa) */
extern  double  u;              /* wind speed (m/sec) */
extern  double  T_g;            /* soil temp at depth z_g (C) */
extern  double  ro;             /* measured runoff (m/sec) */


/*   "sno" update file vars     */

extern  double  time_s;         /* time since start of model run (hrs) */
extern  double  z_s;            /* total snowcover thickness (m) */
extern  double  rho;            /* average snowcover density (kg/m^3) */
extern  double  T_s_0;          /* active snow layer temp (C) */
extern  double  T_s;            /* average snowcover temp (C) */
extern  double  h2o_sat;        /* % of liquid H2O saturation (relative water
				     content, i.e., ratio of water in snowcover
				     to water that snowcover could hold at
				     saturation) */

/*   other snow vars		*/

extern  int     layer_count;    /* number of layers in snowcover: 0, 1, or 2 */
extern  double  z_s_0;          /* active layer depth (m) */
extern  double  z_s_l;          /* lower layer depth (m) */
extern  double  T_s_l;          /* lower layer temp (C) */

extern  double  h2o_vol;        /* liquid h2o content as volume ratio:
				     V_water/(V_snow - V_ice) (unitless) */
extern  double  h2o;            /* liquid h2o content as specific mass
				     (kg/m^2) */
extern  double  h2o_max;        /* max liquid h2o content as specific mass
				     (kg/m^2) */
extern  double  h2o_total;      /* total liquid h2o: includes h2o in snowcover,
				     melt, and rainfall (kg/m^2) */


/*   "mh" update file vars      */

extern  double  time_z;         /* time since start of model run (hrs) */
extern  double  z_u;            /* height of wind measurement (m) */
extern  double  z_T;            /* height of air temp & vapor press meas (m) */
extern  double  z_0;            /* roughness length */


/*   "pre" update file vars     */

extern  double  time_pp;        /* time since start of model run (hrs) */
extern  double  m_pp_data;	/* specific mass of total precip (kg/m^2) */
extern  double  percent_snow;	/* % of total mass that's snow (0 to 1.0) */
extern  double  rho_snow;       /* density of snowfall (kg/m^3) */
extern  double  T_pp;           /* precip temp (C) */

/*   other precip vars 	*/

extern	double	m_rain_data;	/* spec. mass of rain in precip (kg/m^2) */
extern	double	m_snow_data;	/* spec. mass of snow in precip (kg/m^2) */
extern	double	z_snow_data;	/* depth of snow in precip (kg/m^2) */

extern	double	m_pp;		/* m_pp_data adjusted for current timestep */
extern	double	m_rain;		/* m_rain_data   "     "     "       "     */
extern	double	m_snow;		/* m_snow_data   "     "     "       "     */
extern	double	z_snow;		/* z_snow_data   "     "     "       "     */

extern	double	T_rain;		/* rain's temp (K) */
extern	double	T_snow;		/* snowfall's temp (K) */
extern  double  h2o_sat_snow;   /* snowfall's % of liquid H2O saturation */


/*   energy balance vars for current timestep        */

extern  double  R_n;            /* net allwave radiation (W/m^2) */
extern  double  H;              /* sensible heat xfr (W/m^2) */
extern  double  L_v_E;          /* latent heat xfr (W/m^2) */
extern  double  G;              /* heat xfr by conduction & diffusion from soil
				     to snowcover (W/m^2) */
extern  double  G_0;            /* heat xfr by conduction & diffusion from soil
				     or lower layer to active layer (W/m^2) */
extern  double  M;              /* advected heat from precip (W/m^2) */
extern  double  delta_Q;        /* change in snowcover's energy (W/m^2) */
extern  double  delta_Q_0;      /* change in active layer's energy (W/m^2) */

/*   averages of energy balance vars since last output record   */

extern	double	R_n_bar;
extern	double	H_bar;
extern	double	L_v_E_bar;
extern	double	G_bar;
extern	double	G_0_bar;
extern	double	M_bar;
extern	double	delta_Q_bar;
extern	double	delta_Q_0_bar;


/*   mass balance vars for current timestep        */

extern  double  m_s_0;          /* active layer specific mass (kg/m^2) */
extern  double  m_s_l;          /* lower layer specific mass (kg/m^2) */
extern  double  m_s;            /* snowcover's specific mass (kg/m^2) */
extern  double  cc_s;           /* snowcover's cold content (J/m^2) */
extern  double  cc_s_0;         /* active layer cold content (J/m^2) */
extern  double  cc_s_l;         /* lower layer cold content (J/m^2) */
extern  double  melt;       	/* specific melt (kg/m^2 or m) */
extern  double  E;		/* mass flux by evap into air from active
				     layer (kg/m^2/s) */
extern  double  E_s;		/* mass of evap into air & soil from snowcover
				     (kg/m^2) */
extern  double  ro_predict;     /* predicted specific runoff (m/sec) */

/*   sums of some mass balance vars since last output record   */

extern	double	melt_sum;
extern	double	E_s_sum;
extern	double	ro_pred_sum;


/*   work vars                  */

extern  int     rec_count;      /* record counter */

extern  tstep_rec  tstep_info[];	/* array of timestep info:
					   0 : data timestep
					   1 : normal run timestep
					   2 : medium  "     "
					   3 : small   "     "
					 */

extern	double	time_step;	/* current timestep (sec) */

extern  double  start_time;     /* start time (sec) */
extern  double  current_time;   /* time of current time step (sec) */
extern  double  curr_time_hrs;  /* time of current time step (hrs) */

extern	double	time_since_out;	/* time since last output record (sec) */


/*   flags                      */

extern  int     in_file;        /* input records in file? */
extern	int	sn_file;	/* snow properties records in file? */
extern  int     more_sn_recs;   /* any more snow-properties records left? */
extern	int	mh_file;	/* measurement-heights records in file? */
extern  int     more_mh_recs;   /* any more meas-heights records left? */
extern	int	precip_data;	/* precipitation data? */
extern	int	pr_file;	/* precipitation records in file? */
extern  int     more_pr_recs;   /* any more precip records left? */
extern	int	precip_now;	/* precipitation occur for current timestep? */
extern  int     ro_data;        /* runoff data? */
extern  int     out_file;       /* output file specified? */

extern  int     isothermal;     /* melting? */
extern  int     ripe;           /* snow saturated with h2o, ready for runoff? */

extern  int     snowcover;      /* snow on gnd at start of current timestep? */
extern	int	run_no_snow;	/* run model when no snow left? */
extern	int	stop_no_snow;	/* stop model run because no snow left? */
extern	int	temps_in_C;	/* temperatures in degrees C? */

extern	double beta;
extern	double alpha;


/*      declare external subroutines & functions        */
/*      which are part of snobal        */

extern void     EXFUN(adjust_layers,	(NOARGS));
extern void	EXFUN(adjust_snow,	(double delta_z_s, double delta_m_s));
extern void     EXFUN(advec,		(NOARGS));
extern int	EXFUN(below_threshold,	(double threshold));
extern void     EXFUN(calc_layers,	(NOARGS));
extern void	EXFUN(check_fval, 	(float value, float min, float max,
					 char * descrip));
extern double	EXFUN(cold_content,	(double temp, double mass));
extern void	EXFUN(divide_tstep, 	(tstep_rec * tstep));
extern void     EXFUN(do_data_tstep,	(NOARGS));
extern void	EXFUN(do_tstep, 	(tstep_rec * tstep));
extern void     EXFUN(e_bal,		(NOARGS));
extern void     EXFUN(evap_cond,	(NOARGS));
extern int      EXFUN(get_in_rec,	(NOARGS));
extern void     EXFUN(get_mh_rec,	(NOARGS));
extern void     EXFUN(get_pr_rec,	(NOARGS));
extern void     EXFUN(get_sn_rec,	(NOARGS));
extern void     EXFUN(get_args,	  	(int argc, char **argv));
extern void     EXFUN(h2o_compact,	(NOARGS));
extern void     EXFUN(h_le,		(NOARGS));
extern void     EXFUN(init_snow,	(NOARGS));
extern void     EXFUN(initialize,	(NOARGS));
extern void     EXFUN(layer_mass,	(NOARGS));
extern void     EXFUN(mass_bal,		(NOARGS));
extern void     EXFUN(net_rad,		(NOARGS));
extern void     EXFUN(new_density,	(NOARGS));
extern void     EXFUN(open_files,	(NOARGS));
extern void     EXFUN(output,		(NOARGS));
extern void     EXFUN(precip,		(NOARGS));
extern void     EXFUN(runoff,		(NOARGS));
extern void     EXFUN(sn_err,	  	(double new_z_s, double new_rho,
					 double new_T_s_0, double new_T_s,
					 double new_h2o_sat));
extern void     EXFUN(snowmelt,		(NOARGS));
extern void     EXFUN(time_compact,	(NOARGS));
