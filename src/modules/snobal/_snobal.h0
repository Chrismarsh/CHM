/*
** NAME
**      _snobal.h
**
** DESCRIPTION
**      Private include file for the snobal library.
*/

#ifndef _PRIV_SNOBAL_H_
#define _PRIV_SNOBAL_H_

#include "snobal.h"

/* ------------------------------------------------------------------------ */

/*
 *  Private routines in the snobal library.
 */

extern void         _adj_layers	(void);
extern void	      _adj_snow	(double delta_z_s, double delta_m_s);
extern void              _advec	(void);
extern int	   _below_thold	(double threshold);
extern void        _calc_layers	(void);
extern double	  _cold_content	(double temp, double mass);
extern int	  _divide_tstep	(TSTEP_REC * tstep);
extern int	      _do_tstep	(TSTEP_REC * tstep);
extern int               _e_bal	(void);
extern void          _evap_cond (void);
extern void        _h2o_compact (void);
extern int      	  _h_le (void);
extern void         _layer_mass (void);
extern void           _mass_bal (void);
extern void            _net_rad (void);
extern void        _new_density (void);
extern void         _open_files (void);
extern void             _precip (void);
extern void             _runoff (void);
extern void           _snowmelt (void);
extern void       _time_compact (void);

/* ------------------------------------------------------------------------ */

/*
 *  Global variables that are private (internal) to the snobal library.
 */

extern	INPUT_REC input_deltas[];	/* deltas for climate-input parameters
					   over each timestep */

typedef struct {
		double	  m_pp;		/* total precipitation mass (kg/m^2) */
		double	  m_rain;	/* mass of rain in precip (kg/m^2) */
		double	  m_snow;	/*  "   "  snow "     "   (kg/m^2) */
		double	  z_snow;	/* depth of snow in   "   (m) */
	} PRECIP_REC;

extern  PRECIP_REC precip_info[];	/* array of precip info adjusted for
					   each timestep */

extern	int  computed[];		/* array of flags for each timestep;
					   TRUE if computed values for input
					   deltas and precip arrays */

extern	double  h2o_sat_snow;   /* snowfall's % of liquid H2O saturation */
extern	int	isothermal;	/* melting? */
extern	int     snowcover;      /* snow on gnd at start of current timestep? */

#endif /* _PRIV_SNOBAL_H_ */
