//---------------------------------------------------------------------------




#include "_evap_cond.h"
#include "_adj_snow.h"

//---------------------------------------------------------------------------



/*
** NAME
**      _evap_cond -- calculates evaporation/condensation at a point
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_evap_cond(void)
**
** DESCRIPTION
**      Calculates mass lost or gained by evaporation/condensation
**      at a point for 2-layer energy balance snowmelt model snobal.c;
**      Also adjusts the liq h2o, mass and depth of the snow layer;
**      Assumes that liq h2o is favored in evap as the ratio of 
**      vaporization to sublimation (0.882); Half the ice lost as evap
**      is assumed to be lost depth; the rest reduces the density;
**
** GLOBAL VARIABLES READ
**	E	
**	layer_count
**	P_a
**	rho
**	T_g
**	T_s_l
**	T_s_0
**	time_step
**	z_g
**
** GLOBAL VARIABLES MODIFIED
**	E_s h2o_total
**
**	(and those variables modified by "_adj_snow")
*/

#include        "ipw.h"
#include        "_snobal.h"
#include        "envphys.h"

#define VAP_SUB (2.501 / 2.835) /* ratio vaporization to sublimatin */

void
_evap_cond(void)
{
        double  E_s_0;          /* mass of evaporation to air (kg/m^2) */
        double  E_s_l;          /* mass of evaporation to soil (kg/m^2) */
        double  E_l;		/* mass flux by evap/cond to soil (kg/m^2/s) */
        double  e_g;            /* soil vapor press */
        double  e_s_l;          /* lower snow layer's vapor press */
        double  k;              /* soil diffusion coef */
        double  prev_h2o_tot;	/* previous value of h2o_total variable */
        double  q_delta;        /* difference between snow & soil spec hum's */
        double  q_g;            /* soil spec hum */
        double  q_s_l;          /* lower snow layer's spec hum */
        double  rho_air;        /* air density */
        double  T_bar;          /* snow-soil mean temp */

        /*      calculate evaporation or condensation   */

	/*
	 *  If no snow on ground at start of timestep, then just exit.
	 */
	if (!snowcover) {
		E_s = 0.0;
		return;
	}
 
	/*
	 *  Total mass change due to evap/cond at surface during timestep
	 */
        E_s_0 = E * time_step;

	/*
	 *  Adjust total h2o for evaporative losses
	 */
        prev_h2o_tot = h2o_total;

        if (h2o_total > 0.0) {
                h2o_total += (E_s_0 * VAP_SUB);
                if (h2o_total <= 0.0)
                        h2o_total = 0.0;
        }

        /*
	 *  Determine total mass change due to evap/cond at soil
	 */
	if (layer_count == 0) 
		E_s_l = 0.0;
	else {
		if (layer_count == 2) {
       		 	e_s_l = sati(T_s_l);
        		T_bar = (T_g + T_s_l) / 2.0;
		}
		else {  /* layer_count == 1 */
       		 	e_s_l = sati(T_s_0);
       		 	T_bar = (T_g + T_s_0) / 2.0;
		}

		q_s_l = SPEC_HUM(e_s_l, P_a);
		e_g = sati(T_g);
		q_g = SPEC_HUM(e_g, P_a);
		q_delta = q_g - q_s_l;
		rho_air = GAS_DEN(P_a, MOL_AIR, T_bar);
		k = DIFFUS(P_a, T_bar);

		E_l = EVAP(rho_air, k, q_delta, z_g);

	        /* total mass of evap/cond for time step */
		E_s_l = E_l * time_step;

		/** adjust h2o_total for evaporative losses **/
		if (h2o_total > 0.0) {
                	h2o_total += (E_s_l * VAP_SUB);
                	if (h2o_total <= 0.0)
                        	h2o_total = 0.0;
        	}
	}

        E_s = E_s_0 + E_s_l;

        /*      adj mass and depth for evap/cond        */

	if (layer_count > 0)
		_adj_snow( ((E_s + (prev_h2o_tot - h2o_total)) / rho) / 2.0,
			    E_s);
}
