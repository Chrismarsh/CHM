//01/26/07
//---------------------------------------------------------------------------

#ifndef _snobalH
#define _snobalH

    	/*
	 *  default for snowcover's maximum liquid h2o content as volume
	 *  ratio: V_water/(V_snow - V_ice)
	 */
#define DEFAULT_MAX_H2O_VOL  0.01

	/*
	 *  default for maximum active (surface) layer depth (m)
	 */
#define DEFAULT_MAX_Z_S_0  0.25

	/*
	 *  default for depth of soil temperature measurement (m)
	 */
#define DEFAULT_Z_G      0.5

	/*
	 *  Minimum valid snow temperature (C).  This is also what temperatures
	 *  are set to when there's no snow (instead of 0 K).  This yields a
	 *  smaller quantization range in the output image: -75 C to 0 C
	 *  (instead of -273.16 C to 0 C).
	 */
#define MIN_SNOW_TEMP	-75


	/*
	 *  default for medium run timestep (minutes)
	 */
#define	DEFAULT_MEDIUM_TSTEP  15

	/*
	 *  default for small run timestep (minutes)
	 */
#define	DEFAULT_SMALL_TSTEP  1

	/*
	 *  default for normal run timestep's threshold for a layer's mass
	 *  (kg/m^2)
	 */
#define	DEFAULT_NORMAL_THRESHOLD  60.0

	/*
	 *  default for medium run timestep's threshold for a layer's mass
	 *  (kg/m^2)
	 */
#define	DEFAULT_MEDIUM_THRESHOLD  10.0

	/*
	 *  default for small run timestep's threshold for a layer's mass
	 *  (kg/m^2)
	 */
#define	DEFAULT_SMALL_THRESHOLD  1.0

	/*
	 *  Does a time fall within the current input data timestep? 
	 */
#define IN_CURR_DATA_TSTEP(time)	\
		((current_time <= (time)) && \
		 ((time) < current_time + tstep_info[DATA_TSTEP].time_step))

typedef struct {
		double	  m_pp;		/* total precipitation mass (kg/m^2) */
		double	  m_rain;	/* mass of rain in precip (kg/m^2) */
		double	  m_snow;	/*  "   "  snow "     "   (kg/m^2) */
		double	  z_snow;	/* depth of snow in   "   (m) */
	} PRECIP_REC;

extern  PRECIP_REC precip_info[];	/* array of precip info adjusted for each timestep */

                                         
extern	int  computed[];		/* array of flags for each timestep;
					   TRUE if computed values for input
					   deltas and precip arrays */

extern	double  h2o_sat_snow;   /* snowfall's % of liquid H2O saturation */
extern	int	isothermal;	/* melting? */
extern	int     snowcover;      /* snow on gnd at start of current timestep? */

/*   variables that control model execution   */

extern	int	run_no_snow;	/* continue model even if snow disappears? */
extern	int	stop_no_snow;	/* stopped model because no snow left? */

extern	void	(*out_func)(void);	/* -> output function */


/*   constant model parameters  */

extern  double  max_z_s_0;      /* maximum active layer thickness (m) */
extern  double  max_h2o_vol;    /* max liquid h2o content as volume ratio:
				     V_water/(V_snow - V_ice) (unitless) */


/*   time step information */

typedef struct {
		int 	  level;	/* timestep's level */
#define	DATA_TSTEP	   0
#define	NORMAL_TSTEP	   1
#define	MEDIUM_TSTEP	   2
#define	SMALL_TSTEP	   3

		double    time_step;	/* length of timestep (seconds) */
		int	  intervals;	/* # of these timestep that are in
					   the previous-level's timestep
					   (not used for level 0: data tstep) */
		double	  threshold;	/* mass threshold for a layer to use
					   this timestep
					   (not used for level 0: data tstep) */
		int	  output;	/* flags whether or not to call output
					   function for timestep */
#define WHOLE_TSTEP	  0x1		/* output when tstep is not divided */
#define DIVIDED_TSTEP	  0x2		/* output when timestep is divided */

	} TSTEP_REC;

extern  TSTEP_REC  tstep_info[];         /*	 array of info for each timestep:
						   0 : data timestep
						   1 : normal run timestep
						   2 : medium  "     "
						   3 : small   "     "
					 */

extern	double	time_step;	/* length current timestep (sec) */
extern  double  current_time;   /* start time of current time step (sec) */
extern	double	time_since_out;	/* time since last output record (sec) */


/*   snowpack information   */

extern  int     layer_count;    /* number of layers in snowcover: 0, 1, or 2 */
extern  double  z_s;            /* total snowcover thickness (m) */
extern  double  z_s_0;          /* active layer depth (m) */
extern  double  z_s_l;          /* lower layer depth (m) */
extern  double  rho;            /* average snowcover density (kg/m^3) */
extern  double  m_s;            /* snowcover's specific mass (kg/m^2) */
extern  double  m_s_0;          /* active layer specific mass (kg/m^2) */
extern  double  m_s_l;          /* lower layer specific mass (kg/m^2) */
extern  double  T_s;            /* average snowcover temp (K) */
extern  double  T_s_0;          /* active snow layer temp (K) */
extern  double  T_s_l;          /* lower layer temp (C) */
extern  double  cc_s;           /* snowcover's cold content (J/m^2) */
extern  double  cc_s_0;         /* active layer cold content (J/m^2) */
extern  double  cc_s_l;         /* lower layer cold content (J/m^2) */
extern  double  h2o_sat;        /* % of liquid H2O saturation (relative water
				     content, i.e., ratio of water in snowcover
				     to water that snowcover could hold at
				     saturation) */
extern  double  h2o_vol;        /* liquid h2o content as volume ratio:
				     V_water/(V_snow - V_ice) (unitless) */
extern  double  h2o;            /* liquid h2o content as specific mass
				     (kg/m^2) */
extern  double  h2o_max;        /* max liquid h2o content as specific mass
				     (kg/m^2) */
extern  double  h2o_total;      /* total liquid h2o: includes h2o in snowcover,
				     melt, and rainfall (kg/m^2) */


/*   climate-data input records   */

extern  int     ro_data;        /* runoff data? */

typedef struct {
		double S_n;	/* net solar radiation (W/m^2) */
		double I_lw;	/* incoming longwave (thermal) rad (W/m^2) */
		double T_a;	/* air temp (C) */
		double e_a;	/* vapor pressure (Pa) */
		double u;	/* wind speed (m/sec) */
		double T_g;	/* soil temp at depth z_g (C) */
		double ro;	/* measured runoff (m/sec) */
	} INPUT_REC;

extern	INPUT_REC  input_rec1;	/* input data for start of data timestep */
extern	INPUT_REC  input_rec2;	/*   "     "   "  end   "   "      "     */
extern	INPUT_REC input_deltas[];	/* deltas for climate-input parameters over each timestep */

/*   climate-data input values for the current run timestep */

extern  double  S_n;		/* net solar radiation (W/m^2) */
extern  double  I_lw;           /* incoming longwave (thermal) rad (W/m^2) */
extern  double  T_a;            /* air temp (C) */
extern  double  e_a;            /* vapor pressure (Pa) */
extern  double  u;              /* wind speed (m/sec) */
extern  double  T_g;            /* soil temp at depth z_g (C) */
extern  double  ro;             /* measured runoff (m/sec) */


/*   other climate input   */

extern  double  P_a;            /* air pressure (Pa) */


/*   measurement heights/depths   */

extern	int	relative_hts;	/* TRUE if measurements heights, z_T
				   and z_u, are relative to snow
				   surface; FALSE if they are
				   absolute heights above the ground */
extern  double  z_g;            /* depth of soil temp meas (m) */
extern  double  z_u;            /* height of wind measurement (m) */
extern  double  z_T;            /* height of air temp & vapor pressure
				   measurement (m) */
extern  double  z_0;            /* roughness length */


/*   precipitation info for the current DATA timestep    */

extern	int	precip_now;	/* precipitation occur for current timestep? */
extern  double  m_pp;		/* specific mass of total precip (kg/m^2) */
extern  double  percent_snow;	/* % of total mass that's snow (0 to 1.0) */
extern  double  rho_snow;       /* density of snowfall (kg/m^3) */
extern  double  T_pp;           /* precip temp (C) */
extern	double	T_rain;		/* rain's temp (K) */
extern	double	T_snow;		/* snowfall's temp (K) */
extern  double  h2o_sat_snow;   /* snowfall's % of liquid H2O saturation */

/*   precipitation info adjusted for current run timestep   */

extern	double	m_precip;	/* specific mass of total precip (kg/m^2) */
extern	double	m_rain;		/*    "      "   of rain in precip (kg/m^2) */
extern	double	m_snow;		/*    "      "   "  snow "    "    (kg/m^2) */
extern	double	z_snow;		/* depth of snow in precip (m) */


/*   energy balance info for current timestep        */

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

extern  double  melt;       	/* specific melt (kg/m^2 or m) */
extern  double  E;		/* mass flux by evap into air from active
				     layer (kg/m^2/s) */
extern  double  E_s;		/* mass of evap into air & soil from snowcover
				     (kg/m^2) */
extern  double  ro_predict;     /* predicted specific runoff (m/sec) */

/*   sums of mass balance vars since last output record   */

extern	double	melt_sum;
extern	double	E_s_sum;
extern	double	ro_pred_sum;

#endif
 