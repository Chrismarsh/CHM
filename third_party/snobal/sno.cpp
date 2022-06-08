//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

//
// Created by chris on 27/11/15.
//

#include "sno.h"
#include <func/func.hpp>
#include "chm_satw_function.hpp"
#include "chm_sati_function.hpp"
using namespace snobalMacros;

/* FunC TODO:
 * - We can avoid any function evals by reading a good LUT from a json file
 * - higher order interpolation _might_ be desirable */
static FunctionContainer<double> func_w {MyFunction_satw<double>};
static FunctionContainer<double> func_i {MyFunction_sati<double>};

static UniformLinearInterpolationTable<double> satw_lut(&func_w, LookupTableParameters<double> {FREEZE - 50.0, FREEZE + 50.0, 0.00258598});
static UniformLinearInterpolationTable<double> sati_lut(&func_i, LookupTableParameters<double> {90.0, FREEZE, 0.000388715});

//UniformLookupTableGenerator<double> gen_w(&func_w, FREEZE - 50.0, FREEZE + 50.0);
//UniformLookupTableGenerator<double> gen_i(&func_i, 90.0, FREEZE);
//std::string implName = "UniformLinearInterpolationTable";
//// Hard coded stepsizes satisfy a tolerance of 1e-8;
//// gen_by_tol returns unique pointers to LUT:
//static auto satw_lut = gen_w.generate_by_step(implName, 0.00258598);
//static auto sati_lut = gen_i.generate_by_step(implName, 0.000388715);

double sno::ssxfr(
        double k1,    /* layer 1's thermal conductivity (J / (m K sec))  */
        double k2,    /* layer 2's    "         "                        */
        double t1,    /* layer 1's average layer temperature (K)	   */
        double t2,    /* layer 2's    "      "        "         	   */
        double d1,     /* layer 1's thickness (m)			   */
        double d2)     /* layer 2's    "       "			   */
{
    double xfr;

    xfr = 2.0 * (k1 * k2 * (t2 - t1)) / ((k2 * d1) + (k1 * d2));

    return (xfr);
}

double sno::satw(double tk)        /* air temperature (K)		*/
{
     double x;
     double l10;

    if (tk <= 0.)
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("tk < 0"));
    }

    return (satw_lut)(tk);
}

double sno::sati(double tk)        /* air temperature (K)	*/
{
    double l10;
    double x;

    if (tk <= 0.)
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("tk<0"));
    }

    if (tk > FREEZE)
    {
        x = satw(tk);
        return (x);
    }

    if (tk < 90.0)
    {
        return 0.0;
    }

    return (sati_lut)(tk);
}

double sno::new_tsno(
        double spm,    /* layer's specific mass (kg/m^2) 	 */
        double t0,    /* layer's last temperature (K) 	 */
        double ccon)    /* layer's adjusted cold content (J/m^2) */
{
    double tsno;
    double cp;
    double tdif;

    cp = CP_ICE(t0);

    tdif = ccon / (spm * cp);
    tsno = tdif + FREEZE;

    return (tsno);
}

/*
** NAME
        **      init_snow -- initialize the properties for the snowcover
**
** SYNOPSIS
**      #include "snobal.hpp"
**
**      void
        **	init_snow(void)
**
** DESCRIPTION
        **      This routine initializes the properties for the snowcover.  It
        **	determines the number of layers, their individual properties,
        **	the cold content for the snowcover and its layers, and the
**	snowcover's water content.  The following global variables
**	should be initialized before invoking this routine:
**
**		z_s	depth of snowcover (m)
**		rho	density of snowcover (kg/m^3)
**		T_s	average temperature of snowcover (K)
**		T_s_0	temperature of surface layer of snowcover (K)
**		T_s_l	temperature of lower layer of snowcover (K)
**		h2o_sat	% of liquid h2o saturation (0 to 1.0)
**
**		max_h2o_vol	maximum liquid h2o content as volume ratio:
**				    V_water/(V_snow - V_ice) (unitless)
**
** GLOBAL VARIABLES READ
        **	h2o_sat
**	layer_count
        **	m_s_0
**	m_s_l
        **	max_h2o_vol
**	rho
        **	T_s
**	T_s_0
        **	T_s_l
**	z_s
        **
        ** GLOBAL VARIABLES MODIFIED
**	cc_s
        **	cc_s_0
**	cc_s_l
        **	h2o
**	h2o_max
        **	h2o_total
**	h2o_vol
        **	m_s
**	m_s_0
        **	m_s_l
**	rho
        **	T_s
**	T_s_0
        **	T_s_l
*/
void sno::init_snow(void)
{
    double rho_dry;    /* snow density without H2O */

    if (T_s_0 < (MIN_SNOW_TEMP + FREEZE)) // If T_s_0 has gone crazy (very cold), rest to air temperature
    {
        T_s_0 = input_rec1.T_a;
    }

    if (T_s_l < (MIN_SNOW_TEMP + FREEZE)) // If T_s_l has gone crazy (very cold), rest to air temperature
    {
        T_s_l = input_rec1.T_a;
    }

    m_s = rho * z_s;

    _calc_layers();

    if (layer_count == 0)
    {
        /*
         *  If mass > 0, then it must be below threshold.
         *  So turn this little bit of mass into water.
         */
        if (m_s > 0.0)
            h2o_total += m_s;

        rho = 0.0;
        m_s = cc_s = 0.0;
        m_s_0 = cc_s_0 = 0.0;
        m_s_l = cc_s_l = 0.0;

        /*
         *  Note: Snow temperatures are set to MIN_SNOW_TEMP
         *	  (as degrees K) instead of 0 K to keep quantization
         *	  range in output image smaller.
         */
        T_s = T_s_0 = T_s_l = MIN_SNOW_TEMP + FREEZE;
        h2o_vol = h2o = h2o_max = h2o_sat = 0.0;
    }

    else
    {
        /*
         *  Compute specific mass for each layer.
         */
        _layer_mass();

        cc_s_0 = _cold_content(T_s_0, m_s_0);

        if (layer_count == 2)
        {
            cc_s_l = _cold_content(T_s_l, m_s_l);
        }
        else
        {
            T_s_l = MIN_SNOW_TEMP + FREEZE;
            cc_s_l = 0.0;
        }

        /*
         *  Compute liquid water content as volume ratio, and
         *  snow density without water.
         */
        h2o_vol = h2o_sat * max_h2o_vol;
        rho_dry = DRY_SNO_RHO(rho, h2o_vol);

        /*
         *  Determine maximum liquid water content (as specific mass)
         *  and the actual liquid water content (as specific mass).
         */
        h2o_max = H2O_LEFT(z_s, rho_dry, max_h2o_vol);
        h2o = h2o_sat * h2o_max;
    }
}

/* ----------------------------------------------------------------------- */

/*
 * psi-functions
 *	code =	SM	momentum
 *		SH	sensible heat flux
 *		SV	latent heat flux
 */

double sno::psi(
        double zeta,        /* z/lo				*/
        int code)        /* which psi function? (see above) */
{
    double x;        /* height function variable	*/
    double result;

    if (zeta > 0)
    {        /* stable */
        if (zeta > 1)
            zeta = 1;
        result = -BETA_S * zeta;
    }

    else if (zeta < 0)
    {    /* unstable */

        x = sqrt(sqrt(1 - BETA_U * zeta));

        switch (code)
        {
            case SM:
                result = 2 * log((1 + x) / 2) + log((1 + x * x) / 2) -
                         2 * atan(x) + M_PI_2;
                break;

            case SH:
            case SV:
                result = 2 * log((1 + x * x) / 2);
                break;

            default: /* shouldn't reach */
                BOOST_THROW_EXCEPTION(module_error() << errstr_info ("psi-function code not of these: SM, SH, SV"));
                break;
//			bug("psi-function code not of these: SM, SH, SV");
        }
    }

    else
    {            /* neutral */
        result = 0;
    }

    return (result);
}

/* ----------------------------------------------------------------------- */

int sno::hle1(
        double press,    /* air pressure (Pa)			*/
        double ta,    /* air temperature (K) at height za	*/
        double ts,    /* surface temperature (K)		*/
        double za,    /* height of air temp measurement (m)	*/
        double ea,    /* vapor pressure (Pa) at height zq	*/
        double es,    /* vapor pressure (Pa) at surface	*/
        double zq,    /* height of spec hum measurement (m)	*/
        double u,    /* wind speed (m/s) at height zu	*/
        double zu,    /* height of wind speed measurement (m)	*/
        double z0,    /* roughness length (m)			*/

        /* output variables */

        double *h,    /* sens heat flux (+ to surf) (W/m^2)	*/
        double *le,    /* latent heat flux (+ to surf) (W/m^2)	*/
        double *e)    /* mass flux (+ to surf) (kg/m^2/s)	*/
{
    double ah = AH;
    double av = AV;
    double cp = CP_AIR;
    double d0;    /* displacement height (eq. 5.3)	*/
    double dens;    /* air density				*/
    double diff;    /* difference between guesses		*/
    double factor;
    double g = GRAVITY;
    double k = VON_KARMAN;
    double last;    /* last guess at lo			*/
    double lo;    /* Obukhov stability length (eq. 4.25)	*/
    double ltsh;    /* log ((za-d0)/z0)			*/
    double ltsm;    /* log ((zu-d0)/z0)			*/
    double ltsv;    /* log ((zq-d0)/z0)			*/
    double qa;    /* specific humidity at height zq	*/
    double qs;    /* specific humidity at surface		*/
    double ustar;    /* friction velocity (eq. 4.34')	*/
    double xlh;    /* latent heat of vap/subl		*/
    int ier;    /* return error code			*/
    int iter;    /* iteration counter			*/

    /*
     * check for bad input
     */

    /* heights must be positive */
    if (z0 <= 0 || zq <= z0 || zu <= z0 || za <= z0)
    {

        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("height not postive z0=" + std::to_string(z0)
                                                            +" zq="+std::to_string(zq)
                                                             +" zu="+std::to_string(za)

                              ));
//		usrerr ("height not positive; z0=%f\tzq=%f\tzu=%\tza=%f",
//		       z0, zq, zu, za);
        ier = -2;
        return (ier);
    }

    /* temperatures are Kelvin */
    if (ta <= 0 || ts <= 0)
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Temps not K"));
//		usrerr ("temps not K; ta=%f\tts=%f", ta, ts);
        ier = -2;
        return (ier);
    }

    /* pressures must be positive */
    if (ea <= 0 || es <= 0 || press <= 0 || ea >= press || es >= press)
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Press <0"));
//		usrerr ("press < 0; ea=%f\tes=%f\tpress=%f", ea, es, press);
        ier = -2;
        return (ier);
    }

    /* vapor pressures can't exceed saturation */
    /* if way off stop */
    if ((es - 25.0) > sati(ts) || (ea - 25.0) > satw(ta))
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("vp > sat"));
//		usrerr ("vp > sat; es=%f\tessat=%f\tea=%f\teasat=%f",
//			es, sati(ts), ea, sati(ta));
        ier = -2;
        return (ier);
    }
    /* else fix them up */
    if (es > sati(ts))
    {
        es = sati(ts);
    }
    if (ea > satw(ta))
    {
        ea = satw(ta);
    }

    /*
     * displacement plane height, eq. 5.3 & 5.4
     */

    d0 = 2 * PAESCHKE * z0 / 3;

    /*
     * constant log expressions
     */

    ltsm = log((zu - d0) / z0);
    ltsh = log((za - d0) / z0);
    ltsv = log((zq - d0) / z0);

    /*
     * convert vapor pressures to specific humidities
     */
    qa = SPEC_HUM(ea, press);
    qs = SPEC_HUM(es, press);

    /*
     * convert temperature to potential temperature
     */

    ta += DALR * za;

    /*
     * air density at press, virtual temp of geometric mean
     * of air and surface
     */

    dens = GAS_DEN(press, MOL_AIR,
                   VIR_TEMP(sqrt(ta * ts), sqrt(ea * es), press));

    /*
     * starting value, assume neutral stability, so psi-functions
     * are all zero
     */

    ustar = k * u / ltsm;
    factor = k * ustar * dens;
    *e = (qa - qs) * factor * av / ltsv;
    *h = (ta - ts) * factor * cp * ah / ltsh;

    /*
     * if not neutral stability, iterate on Obukhov stability
     * length to find solution
     */

    iter = 0;
    if (ta != ts)
    {

        lo = HUGE_VAL;

        do
        {
            last = lo;

            /*
             * Eq 4.25, but no minus sign as we define
             * positive H as toward surface
             */

            /*
             * There was an error in the old version of this
             * line that omitted the cubic power of ustar.
             * Now, this error has been fixed.
             */

            lo = ustar * ustar * ustar * dens
                 / (k * g * (*h / (ta * cp) + 0.61 * *e));

            /*
             * friction velocity, eq. 4.34'
             */

            ustar = k * u / (ltsm - psi(zu / lo, SM));

            /*
             * evaporative flux, eq. 4.33'
             */

            factor = k * ustar * dens;
            *e = (qa - qs) * factor * av /
                 (ltsv - psi(zq / lo, SV));

            /*
             * sensible heat flus, eq. 4.35'
             * with sign reversed
             */

            *h = (ta - ts) * factor * ah * cp /
                 (ltsh - psi(za / lo, SH));

            diff = last - lo;

        } while (fabs(diff) > THRESH &&
                 fabs(diff / lo) > THRESH &&
                 ++iter < ITMAX);
    }

    if(iter >= ITMAX || std::isinf(diff) ) // fall back to neutral
    {
        ier = -1;

        //failed to converge, likely low winds, assume neutral
        ustar = k * u / ltsm;
        factor = k * ustar * dens;
        *e = (qa - qs) * factor * av / ltsv;
        *h = (ta - ts) * factor * cp * ah / ltsh;
        ier = 0;
    }
    else
    {
        ier = 0;
    }



//    ier = (iter >= ITMAX) ? -1 : 0;

    xlh = LH_VAP(ts);
    if (ts <= FREEZE)
        xlh += LH_FUS(ts);

    /*
     * latent heat flux (- away from surf)
     */
    *le = xlh * *e;

    return (ier);
}

double sno::heat_stor(
        double cp,    /* specific heat of layer (J/kg K) */
        double spm,    /* layer specific mass (kg/m^2)    */
        double tdif)    /* temperature change (K)          */
{
    double stor;

    stor = cp * spm * tdif;

    return (stor);
}

double sno::g_soil(
        double rho,    /* snow layer's density (kg/m^3)	     */
        double tsno,    /* snow layer's temperature (K)		     */
        double tg,    /* soil temperature (K)			     */
        double ds,    /* snow layer's thickness (m)		     */
        double dg,    /* dpeth of soil temperature measurement (m) */
        double pa)    /* air pressure (Pa)			     */
{
    double k_g;
    double kcs;
    double k_s;
    double g;

/*	check tsno	*/
    if (tsno > FREEZE)
    {
//		warn("g_soil: tsno = %8.2f; set to %8.2f\n", tsno, FREEZE);
        tsno = FREEZE;
    }

/*	set effective soil conductivity	*/
    k_g = efcon(KT_WETSAND, tg, pa);

/*	calculate G	*/
    /*	set snow conductivity	*/
    kcs = KTS(rho);
    k_s = efcon(kcs, tsno, pa);

    g = ssxfr(k_s, k_g, tsno, tg, ds, dg);

    return (g);
}

double sno::g_snow(
        double rho1,    /* upper snow layer's density (kg/m^3)	*/
        double rho2,    /* lower  "     "        "    (kg/m^3)	*/
        double ts1,    /* upper snow layer's temperature (K)	*/
        double ts2,    /* lower  "     "         "       (K)	*/
        double ds1,    /* upper snow layer's thickness (m)	*/
        double ds2,    /* lower  "     "         "     (m)	*/
        double pa)    /* air pressure (Pa)			*/
{
    double kcs1;
    double kcs2;
    double k_s1;
    double k_s2;
    double g;


/*	calculate G	*/
    if (ts1 == ts2)
        g = 0.0;
    else
    {
        /*	set snow conductivity	*/
        kcs1 = KTS(rho1);
        kcs2 = KTS(rho2);
        k_s1 = efcon(kcs1, ts1, pa);
        k_s2 = efcon(kcs2, ts2, pa);

        /*	calculate g	*/
        g = ssxfr(k_s1, k_s2, ts1, ts2, ds1, ds2);
    }

    return (g);
}

double sno::efcon(
        double k,    /* layer thermal conductivity (J/(m K sec)) */
        double t,    /* layer temperature (K)		    */
        double p)    /* air pressure (Pa)  			    */
{
    double etc;
    double de;
    double lh;
    double e;
    double q;

    /*	calculate effective layer diffusion
        (see Anderson, 1976, pg. 32)		*/
    de = DIFFUS(p, t);

    /*	set latent heat from layer temp.	*/
    if (t > FREEZE)
        lh = LH_VAP(t);
    else if (t == FREEZE)
        lh = (LH_VAP(t) + LH_SUB(t)) / 2.0;
    else
        lh = LH_SUB(t);

    /*	set mixing ratio from layer temp.	*/
    e = sati(t);
    q = MIX_RATIO(e, p);

    /*	calculate effective layer conductivity	*/
    etc = k + (lh * de * q);

    return (etc);
}

/*
** NAME
**      do_data_tstep -- run model for 1 data timestep between 2 input records
**
** SYNOPSIS
**
**	int
**	do_data_tstep(void)
**
** DESCRIPTION
**	This routine performs the model's calculations for 1 data timestep
**	between 2 input-data records which are in 'input_rec1' and
**	'input_rec2'.
**
**	If there's precipitation during the data timestep, the flag
**	'precip_now' used be TRUE.  Furthermore, the routine requires
**	that the following precipitation variables have been initialized:
**
**		m_pp
**		percent_snow
**		rho_snow
**		T_pp
**
**	This routine divides the data timestep into the appropriate number
**	of normal run timesteps.  The input values for each normal timestep
**	are computed from the two input records by linear interpolation.
**
**	If output is desired for any of the run timesteps (normal, medium,
**	or small), the appropriate output flags must be set in the proper
**	timestep's record (i.e., the array 'tstep_info').  If any output
**	flag is set, the routine requires that the global variable 'out_func'
**	point to appropriate output function.
**
**	This routine may return in the middle of a data timestep if:
**
**		a)  the output function pointed to by 'out_func' is called, and
**		b)  the flag 'run_no_snow' is FALSE, and
**		c)  there is no snow remaining on the ground at the end of
**		    timestep
**
**	In this happens, the flag 'stop_no_snow' is set to TRUE.
**
** RETURN VALUE
**
**	TRUE	The model's calculations were completed.
**
**	FALSE	An error occured, and a message explaining the error has
**		been stored with the 'usrerr' routine.
**
** GLOBAL VARIABLES READ
**	e_a
**	I_lw
**	in_rec
**	layer_count
**	m_pp_data
**	m_rain_data
**	m_snow_data
**	more_pr_recs
**	precip_data
**	ro
**	ro_data
**	run_no_snow
**	S_n
**	T_a
**	T_g
**	tstep_info
**	u
**	z_snow_data
**
** GLOBAL VARIABLES MODIFIED
**	precip_now
**	stop_no_snow
*/

    int sno::do_data_tstep(void)
    {
        PRECIP_REC *pp_info = precip_info; /* precip info for data timestep */
        TSTEP_REC *data_tstep = tstep_info; /* timestep info for data timestep */
        int level;            /* loop index */

/*
*  Copy values from first input record into global variables.
*/
        S_n = input_rec1.S_n;
        I_lw = input_rec1.I_lw;
        T_a = input_rec1.T_a;
        e_a = input_rec1.e_a;
        u = input_rec1.u;
        T_g = input_rec1.T_g;
        if (ro_data)
            ro = input_rec1.ro;

        /*
         *  Compute deltas for the climate input parameters over
         *  the data timestep.
         */
        input_deltas[DATA_TSTEP].S_n = input_rec2.S_n - input_rec1.S_n;
        input_deltas[DATA_TSTEP].I_lw = input_rec2.I_lw - input_rec1.I_lw;
        input_deltas[DATA_TSTEP].T_a = input_rec2.T_a - input_rec1.T_a;
        input_deltas[DATA_TSTEP].e_a = input_rec2.e_a - input_rec1.e_a;
        input_deltas[DATA_TSTEP].u = input_rec2.u - input_rec1.u;
        input_deltas[DATA_TSTEP].T_g = input_rec2.T_g - input_rec1.T_g;
        if (ro_data)
            input_deltas[DATA_TSTEP].ro = input_rec2.ro - input_rec1.ro;

        /*
         *  If there is precipitation, then compute the amount of rain &
         *  snow in it.
         */
        if (precip_now)
        {
            pp_info->m_pp = m_pp;
            pp_info->m_snow = percent_snow * m_pp;
            pp_info->m_rain = m_pp - pp_info->m_snow;
            if (pp_info->m_snow > 0.0)
            {
                if (rho_snow > 0.0)
                    pp_info->z_snow = pp_info->m_snow / rho_snow;
                else
                {
                    LOG_ERROR <<  "rho_snow is <= 0.0 with %_snow > 0.0";
                    return 0;
                }
            }
            else
                pp_info->z_snow = 0.0;

            /*
             *  Mixed snow and rain
             */
            if ((pp_info->m_snow > 0.0) && (pp_info->m_rain > 0.0))
            {
                T_snow = FREEZE;
                h2o_sat_snow = 1.0;
                T_rain = T_pp;
            }

                /*
                 *  Snow only
                 */
            else if (pp_info->m_snow > 0.0)
            {
                if (T_pp < FREEZE)
                {        /* Cold snow */
                    T_snow = T_pp;
                    h2o_sat_snow = 0.0;
                }
                else
                {                /* Warm snow */
                    T_snow = FREEZE;
                    h2o_sat_snow = 1.0;
                }
            }

                /*
                 *  Rain only
                 */
            else if (pp_info->m_rain > 0.0)
            {
                T_rain = T_pp;
            }
        }

        /*
         *  Clear the 'computed' flag at the other timestep levels.
         */
        for (level = NORMAL_TSTEP; level <= SMALL_TSTEP; level++)
            computed[level] = 0;

        /*
         *  Divide the data timestep into normal run timesteps.
         */
        return _divide_tstep(data_tstep);
    }

/*
** NAME
**      _time_compact_ori -- compact snowcover by gravity over time (original param.)
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_time_compact_ori(void)
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
void sno::_time_compact_ori(void)
{
    /*
     *  Maximum density due to compaction by gravity (kg/m^2).
     */
    int A = 350;
    /*
     *  Time when half "saturation", i.e., maximum density is reached
     *  (seconds).
     *  (864000 = 10 days * 24 hours/day * 60 mins/hr * 60 secs/min)
     */

    int B = 86400;
    double time;    /* point on time axis corresponding to current
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
    rho = A / (1 + B / (time + time_step));

    /*
 *  Adjust the snowcover for this new density.
 */
    _new_density();
}

/* 
** NAME 
**      _time_compact -- snowcover gravety, depth & temperature compaction 
** 
** SYNOPSIS 
**      #include "_snobal.h" 
** 
**      void 
**	_time_compact(void) 
** 
** DESCRIPTION 
**	This routine replaces the original simple gravety compaction routine 
**	which "aged" the snowcover by accounting for the compaction
**	or densification by gravity as time passes.  The original routine 
**	relied on an approximation based on equations in Anderson (1976), which 
**	increased the snowcover's density using the following "half-saturation" 
**	function: 
** 
**		rho(time) = A / (1 + B/time) 
** 
**	Historically, precise snow density was not a concern, as long as mass 
**	and SWE were correct.  With the development of the ASO program providing 
**	time-series lidar snow depth images, and the 2017 snow season in the 
*8	southern Sierra Nevada, with individual storms approaching 500mm of 
**	deposition in a few days, and upper elevation snow depths of greater 
**	than 10m, it became clear that a more robust density model was required. 
** 
**	Snow Density:  Snow density is initially a function of the temperature 
**	of the ice particles (snow flakes) as they fall to the ground during a 
**	storm.  Under very cold conditions (-10 to -15 C), we can get snow as 
**	light as 50 kg m-3, which is really light powder and great skiing. 
**	As the ice particle temperature approaches 0C, new snow densities can 
**	be as high as 200 kg m-3, which is not light powder, but still pretty 
**	good skiing, unless there is a lot of it. There are several (4) 
**	processes that can increase snow density during and after the storm. 
**	Note that the largest and most rapid changes occur during or just 
**	following the storm.  Compaction (an increase in density without a 
**	change in mass) can be caused by: 
**	
**	1) Destructive mechanical metamorphism (compaction due to wind - 
**	   affects mainly low density near-surface or new snow) 
**	2) Destructive temperature metamorphism 
**	3) Pressure - compaction due to snow load or overburden (affects both 
**	   new snow deposition, snow on the ground, including drifting) 
**	4) Addition of liquid water due to melt or rain 
** 
**	Of these compaction processes, iSnobal accounts three - 2,3 & 4. 
**	This routine addresses 2, temperature metamorphism, and 3., overburden 
**	metamorphism. The 4th, addition of liquid water is accounted for in the 
**	routine _h2o_compact.c. We are using equations found in Anderson (1976) 
**	and Oleson, et al. (2013), who got his equations from Anderson in the 
**	first place.  (Oleson made it easier to figure out the units...) 
** 
**	Though many dedicated individuals have worked on the density issue over 
**	over the last 60+ years, we have, in general, a group working in Saporro 
**	Japan to thank for most of what we know about snow density. Anderson 
**	got the data and basic fitting equations from careful field and cold 
**	room measurements made by Yosida (1963), Mellor (1964) and Kojima (1967). 
**	The equations we are using are based on those that they derived from data 
**	that they carefully collected over a number of years.  It is noteworthy 
**	that while rapidly changing computers have made the kind of spatial 
**	modeling we are attempting possible, snow physics remains unchanged – and 
**	the field and labortory efforts and data fitting equations from more than 
**	half a century ago represent our best understanding of those physics. 
** 
**	Tz = 0.0 (freezing temperature, C or K) 
**	Ts = snow or precipitation temperature (C or K) 
**	rho = intital snow density (kg/(m^3)) 
**	SWE = snow mass (mm/(m^2)) 
**	K = temperature metamorphism coef. 
**	rho_n = new snow density (kg/(m^3)) 
**	zs = new snow depth (mm) 
** 
**	Proportional Destructive Temperature Metamorphism (PTM): 
** 
**	if (rho < 100) 
**		K = 1.0 
**	else 
**		K = exp(-0.046 * (rho - 100)) 
**	
**	PTM = 0.01 * K * exp(-0.04 * (Tz - Ts)) 
** 
**	Proportional Overburden Compaction (POC): 
** 
**	POC = (0.026 * exp(-0.08 * (Tz - Ts)) * SWE * exp(-21.0 * rho)) 
** 
**	New snow density and depth 
**	
**	rho_n = rho + ((PTM + POC) * rho) 
**	zs_n = SWE / rho_n 
**      
** GLOBAL VARIABLES READ 
**	time_step, T_s, m_s, rho 
** 
** GLOBAL VARIABLES MODIFIED 
**	rho 
** 
*/

void sno::_time_compact(void) 
{
	double	c11;	/* temperature metamorphism coefficient (Anderson, 1976) */
        double	Tz;	/* Freezing temperature (K) */ 	
        double	d_rho_m;
        double	d_rho_c;
        double	rate;

 	/*
 	 *  If the snow is already at or above the maximum density due to
 	 *  compaction, then just leave.
 	 */
 	if ((!snowcover) || (rho >= RMX))
 		return;

 	Tz = FREEZE;
 	/*
 	 *  Calculate rate which compaction will be applied per time step.
 	 *  Rate will be adjusted as time step varies.
 	 */
 	if (m_s >= SWE_MAX)
 		rate = 1.0;
 	else {
 		rate = RD1 * cos((PI * m_s) / SWE_MAX) + RD2;
 		rate = rate / (time_step / nsec_hour);
 	} 

	/** Proportional Destructive Temperature Metamorphism (d_rho_m) **/

 	if (rho < 100)
 		c11 = 1.0;
 	else
 		c11 = exp(-0.046 * (rho - 100));

 	d_rho_m = 0.01 * c11 * exp(-0.04 * (Tz - T_s));
 	d_rho_m /= rate;

 	/** Proportional Overburden Compaction (d_rho_c) **/

 	double slope_ms = m_s;
 	if(slope > 0) // if we have a slope here, apply a cosine correction for compaction
            slope_ms /= std::max(0.001,cos(slope));

 	d_rho_c = (0.026 * exp(-0.08 * (Tz - T_s)) * slope_ms * exp(-21.0 * (rho / water)));
 	d_rho_c /= rate;

	/**	Compute New snow density	**/

 	rho = rho + ((d_rho_m + d_rho_c) * rho);

         /*
 	 *  Adjust the snowcover for this new density.
 	 */
 	_new_density();
}


/*
** NAME
**      _snowmelt -- calculates melt or re-freezing at a point
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_snowmelt(void)
**
** DESCRIPTION
**      Calculates melting or re-freezing for point 2-layer energy balance
**      snowmelt model.
**
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**
*/



    void  sno::_snowmelt(void)
    {
        double Q_0;            /* energy available for surface melt */
        double Q_l;        /* energy available for lower layer melt */
        double Q_freeze;       /* energy used for re-freezing */
        double Q_left;         /* energy left after re_freezing */
        double h2o_refrozen;   /* amount of liquid H2O that was refrozen */


        /*
         *  If no snow on ground at start of timestep, then just exit.
         */
        if (!snowcover)
        {
            melt = 0.0;
            return;
        }

        /***    calculate melt or freezing, and adjust cold content */

        /*** calculate surface melt ***/

        /* energy for surface melt */
        Q_0 = (delta_Q_0 * time_step) + cc_s_0;

        if (Q_0 > 0.0)
        {
            melt = MELT(Q_0);
            cc_s_0 = 0.0;
        }
        else if (Q_0 == 0.0)
        {
            melt = 0.0;
            cc_s_0 = 0.0;
        }
        else
        {
            melt = 0.0;
            cc_s_0 = Q_0;
        }


        /*** calculate lower layer melt ***/

        if (layer_count == 2)
        {
            /* energy for layer melt */
            Q_l = ((G - G_0) * time_step) + cc_s_l;

            if (Q_l > 0.0)
            {
                melt += MELT(Q_l);
                cc_s_l = 0.0;
            }
            else if (Q_l == 0.0)
                cc_s_l = 0.0;
            else
                cc_s_l = Q_l;
        }
        else
        {  /* layer_count == 1 */
            Q_l = 0.0;
        }

        h2o_total += melt;


        /*** adjust layers for re-freezing ***/

        /*    adjust surface layer    */

        h2o_refrozen = 0.0;

        if (cc_s_0 < 0.0)
        {
            /* if liquid h2o present, calc refreezing and adj cc_s_0 */
            if (h2o_total > 0.0)
            {
                Q_freeze = h2o_total * (z_s_0 / z_s) * LH_FUS(FREEZE);
                Q_left = Q_0 + Q_freeze;

                if (Q_left <= 0.0)
                {
                    h2o_refrozen = h2o_total * (z_s_0 / z_s);
                    cc_s_0 = Q_left;
                }
                else
                {
                    h2o_refrozen = (h2o_total * (z_s_0 / z_s)) -
                                   MELT(Q_left);
                    cc_s_0 = 0.0;
                }
            }
        }

        /*    adjust lower layer for re-freezing */

        if ((layer_count == 2) && (cc_s_l < 0.0))
        {
            /* if liquid h2o, calc re-freezing and adj cc_s_l */
            if (h2o_total > 0.0)
            {
                Q_freeze = h2o_total * (z_s_l / z_s) * LH_FUS(FREEZE);
                Q_left = Q_l + Q_freeze;

                if (Q_left <= 0.0)
                {
                    h2o_refrozen += h2o_total * (z_s_l / z_s);
                    cc_s_l = Q_left;
                }
                else
                {
                    h2o_refrozen += ((h2o_total * (z_s_l / z_s)) -
                                     MELT(Q_left));
                    cc_s_l = 0.0;
                }
            }
        }

        /*
         *  Note:  because of rounding errors, h2o_refrozen may not
         * 	   be exactly the same as h2o_total.  Check for this
         *	   case, and if so, then just zero out h2o_total.
         */
        if (fabs(h2o_total - h2o_refrozen) <= 1e-8)
        {
            h2o_total = 0.0;
        } else
        {
            h2o_total -= h2o_refrozen;
        }

        /***	determine if snowcover is isothermal    ***/

        double eps = 1e-6;
//        if ((layer_count == 2) && (cc_s_0 == 0.0) && (cc_s_l == 0.0))
        if ((layer_count == 2) && std::fabs(cc_s_0 - 0.0) < eps && std::fabs(cc_s_l - 0.0) < eps)
            isothermal = 1;
        else if ((layer_count == 1) && std::fabs(cc_s_0 - 0.0) < eps)
            isothermal = 1;
        else
            isothermal = 0;

        /***    adjust depth and density for melt  ***/

        if (melt > 0.0)
            _adj_snow(-(melt / rho), 0.0);

        /***    set total cold content   ***/
        if (layer_count == 2)
            cc_s = cc_s_0 + cc_s_l;
        else if (layer_count == 1)
            cc_s = cc_s_0;
    }

/*
** NAME
**      _runoff -- calculates runoff from snowcover
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_runoff(void)
**
** DESCRIPTION
**      Calculates runoff for point energy budget 2-layer snowmelt model
**
** GLOBAL VARIABLES READ
**	h2o_total
**	layer_count
**	snowcover
**	max_h2o_vol
**	z_s
**
** GLOBAL VARIABLES MODIFIED
**	h2o
**	h2o_max
**	h2o_sat
**	h2o_vol
**	rho
**	ro_predict
*/

    void sno::_runoff(void)
    {
        double m_s_dry;    /* snowcover's mass without liquid H2O */
        double rho_dry;    /* snow density without liquid H2O */

        /* calculate runoff */

        /*
         *  If no snow on ground at start of timestep or no layers currently,
         *  then all water (e.g., rain) is runoff.
         */
        if ((!snowcover) || (layer_count == 0))
        {
            ro_predict = h2o_total;
            return;
        }

        /*
     *  Determine the snow density without any water, and the maximum
     *  liquid water the snow can hold.
     */
        m_s_dry = m_s - h2o_total;
        rho_dry = m_s_dry / z_s;
        h2o_max = H2O_LEFT(z_s, rho_dry, max_h2o_vol);

        /*
         *  Determine runoff, and water left in the snow
         */
        if (h2o_total > h2o_max)
        {
            ro_predict = h2o_total - h2o_max;
            h2o = h2o_max;
            h2o_sat = 1.0;
            h2o_vol = max_h2o_vol;

            /*
             *  Update the snowcover's mass for the loss of runoff.
             */
            _adj_snow(0.0, -ro_predict);
        }
        else
        {
            ro_predict = 0.0;
            h2o = h2o_total;
            h2o_sat = h2o / h2o_max;
            h2o_vol = h2o_sat * max_h2o_vol;
        }
    }

/*
** NAME
**      _precip -- process a precipitation event
**
** SYNOPSIS
**	#include "_snobal.h"
**
**      void
**	_precip(void)
**
** DESCRIPTION
**      This routine processes a precipitation event, i.e., the current
**	precip record, if there's one for the current timestep.  It
**	determines if the precip is rain or snow which increases the
**	snowcover.
**
** GLOBAL VARIABLES READ
**	h2o_sat_snow
**	m_rain
**	m_precip
**	max_h2o_vol
**	precip_now
**	rho_snow
**	snowcover
**	T_snow
**	z_snow
**
** GLOBAL VARIABLES MODIFIED
**	h2o
**	h2o_sat
**	h2o_total
**	rho
**	T_s
**	T_s_0
**	z_s
*/

void sno::_precip(void)
{
    double h2o_vol_snow;    /* liquid water content of new snowfall as
				   volume ratio */

    if (precip_now)
    {
        if (snowcover)
        {
            /*
             *  Adjust snowcover's depth and mass by snowfall's
             *  depth and the total precipitation mass.
             */
            _adj_snow(z_snow, m_precip);

            /*
             *  Determine the additional liquid water that's in
             *  the snowfall, and then add its mass to liquid
             *  water in the whole snowcover.
             */
            h2o_vol_snow = h2o_sat_snow * max_h2o_vol;
            h2o += H2O_LEFT(z_snow, rho_snow, h2o_vol_snow);
        }
        else
        {
            /*
             *  Use snowfall, if any, to setup a new snowcover.
             */
            if (m_snow > 0.0)
            {
                z_s = z_snow;
                rho = rho_snow;
                T_s = T_snow;
                T_s_0 = T_snow;
                T_s_l = T_snow;
                h2o_sat = h2o_sat_snow;

                init_snow();
            }
        }

        /*
         *  Add rainfall and water in the snowcover to total
         *  liquid water.
         */
        h2o_total += h2o + m_rain;
    }
    else
    {
        /*
         *  Add water in the snowcover to total liquid water.
         */
        h2o_total += h2o;
    }
}

/*
** NAME
**      _new_density -- adjust the snowcover's depth and layers for new density
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_new_density(void)
**
** DESCRIPTION
**      This routine adjusts the snowcover's depth for a new density.  The
**	layers are also adjusted accordingly.
**
** GLOBAL VARIABLES READ
**	m_s
**	rho
**
** GLOBAL VARIABLES MODIFIED
**	z_s
**
**	(and those variables modified by "_adj_layers")
*/
void sno::_new_density(void)
{
    z_s = m_s / rho;

    _adj_layers();
}


/*
** NAME
**      _net_rad -- calculates net allwave radiation
**
** SYNOPSIS
**
**	#include "_snobal.h"
**
**      void
**	_net_rad(void)
**
** DESCRIPTION
**      Calculates net allwave radiation from the net solar radiation
**	incoming thermal/longwave radiation, and the snow surface
**	temperature.
**
** GLOBAL VARIABLES READ
**	I_lw
**	S_n
**	T_s_0
**
** GLOBAL VARIABLES MODIFIED
**	R_n
*/
void sno::_net_rad(void)
{
    R_n = S_n + (SNOW_EMISSIVITY * (I_lw - STEF_BOLTZ * pow(T_s_0, 4)));
}

/*
** NAME
**      _mass_bal -- calculates point mass budget of 2-layer snowcover
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_mass_bal(void)
**
** DESCRIPTION
**      Calculates the point mass budget for 2-layer energy budget snowmelt
**	model.  It then solves for new snow temperatures.
**
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**
*/
void sno::_mass_bal(void)
{
    /***    adjust mass and calc. runoff    ***/

    /*	age snow by compacting snow due to time passing */
    if(param_snow_compaction == 1)
    {
         _time_compact();
    }
    else if(param_snow_compaction ==0)
    {
         _time_compact_ori();
    }
    

    /*	process precipitation event */
    _precip();

    /*      calculate melt or freezing and adjust cold content */

    _snowmelt();

    /*      calculate evaporation and adjust snowpack       */

    _evap_cond();

    /*	compact snow due to H2O generated (melt & rain) */
    _h2o_compact();

    /*      calculate runoff, and adjust snowcover */

    _runoff();

    /*
 *  adjust layer temps if there was a snowcover at start of the
 *  timestep and there's still snow on the ground
 */
    if (snowcover)
    {
        if (layer_count == 1)
        {
            T_s_0 = new_tsno(m_s_0, T_s_0, cc_s_0);
            T_s = T_s_0;
        }
        else if (layer_count == 2)
        {
            if (isothermal)
                T_s = T_s_l = T_s_0 = FREEZE;
            else
            {
                T_s_0 = new_tsno(m_s_0, T_s_0, cc_s_0);
                T_s_l = new_tsno(m_s_l, T_s_l, cc_s_l);
                T_s = new_tsno(m_s, T_s, cc_s);
            }
        }
    }
}

/*
** NAME
**      _layer_mass -- calculate the specific mass for each snow layer
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_layer_mass(void)
**
** DESCRIPTION
**      This routine computes the specific mass for each snow layer in
**	the snowcover.  A layer's mass is based its depth and the
**	average snowcover density.
**
** GLOBAL VARIABLES READ
**	layer_count
**	rho
**	z_s_0
**	z_s_l
**
** GLOBAL VARIABLES MODIFIED
**	m_s_0
**	m_s_l
*/
void sno::_layer_mass(void)
{
    if (layer_count == 0)
    {
        m_s_0 = 0.0;
        m_s_l = 0.0;
    }
    else
    {  /* layer_count is 1 or 2 */
        m_s_0 = rho * z_s_0;
        if (layer_count == 2)
            m_s_l = rho * z_s_l;
        else
            m_s_l = 0.0;
    }
}

/*
** NAME
**      _h_le -- calculates turbulent transfer at a point
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      int
**	_h_le(void)
**
** DESCRIPTION
**      Calculates point turbulent transfer (H and L_v_E) for a 2-layer
**	snowcover.
**
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**
*/
int sno::_h_le(void)
{
    double e_s;
    double sat_vp;
    double rel_z_T;  /* relative z_T (temperature measurement
			     height) above snow surface */
    double rel_z_u;  /* relative z_u (windspeed measurement
			     height) above snow surface */

    /* calculate saturation vapor pressure */

    e_s = sati(T_s_0);

    /*** error check for bad vapor pressures ***/

    sat_vp = sati(T_a);
    if (e_a > sat_vp)
    {
        e_a = sat_vp;
    }

    /* determine relative measurement heights */
    if (relative_hts)
    {
        rel_z_T = z_T;
        rel_z_u = z_u;
    } else
    {
        rel_z_T = z_T - z_s;
        rel_z_u = z_u - z_s;
    }

    /* calculate H & L_v_E */

    if (hle1(P_a, T_a, T_s_0, rel_z_T, e_a, e_s, rel_z_T, u,
             rel_z_u, z_0, &H, &L_v_E, &E) != 0)
    {
        LOG_DEBUG << "hle1 did not converge";// sprintf("hle1 did not converge\nP_a %f, T_a %f, T_s_0 %f\nrelative z_T %f, e_a %f, e_s %f\nu %f, relative z_u %f, z_0 %f\n", P_a, T_a, T_s_0, rel_z_T, e_a, e_s, u, rel_z_u, z_0);
        return 0;
    }

    return 1;
}

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

void sno::_h2o_compact(void)
{
    double MAX_DENSITY   = 550;
    /*
     *  Maximum density due to compaction by liquid H2O added (kg/m^2).
     */

    /* double   B    = 0.05; Old value used in snobal*/
      double   B    = 0.4; /* New value from  https://gitlab.com/ars-snow/ipw/tree/master/src/lib/libsnobal/snobal
                            Provide more realistic evolution of snow density during rain-on snow events. 
                           */  

    /*
     *  ratio where half the difference between maximum density and
     *  current density is reached (ratio from 0.0 to 1.0).
     */

    double A;        /* difference between maximum & current
				   densities */
    double h2o_added;    /* ratio of mass of liquid H2O added from
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
    if (h2o_added > 0.000001)
    {
        rho += A / (1 + B / h2o_added);

        /*
     *  Adjust the snowcover for this new density.
     */
        _new_density();
    }
}


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
void sno::_evap_cond(void)
{
#define VAP_SUB (2.501 / 2.835) /* ratio vaporization to sublimatin */

    double E_s_0;          /* mass of evaporation to air (kg/m^2) */
    double E_s_l;          /* mass of evaporation to soil (kg/m^2) */
    double E_l;        /* mass flux by evap/cond to soil (kg/m^2/s) */
    double e_g;            /* soil vapor press */
    double e_s_l;          /* lower snow layer's vapor press */
    double k;              /* soil diffusion coef */
    double prev_h2o_tot;    /* previous value of h2o_total variable */
    double q_delta;        /* difference between snow & soil spec hum's */
    double q_g;            /* soil spec hum */
    double q_s_l;          /* lower snow layer's spec hum */
    double rho_air;        /* air density */
    double T_bar;          /* snow-soil mean temp */

    /*      calculate evaporation or condensation   */

    /*
     *  If no snow on ground at start of timestep, then just exit.
     */
    if (!snowcover)
    {
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

    if (h2o_total > 0.0)
    {
        h2o_total += (E_s_0 * VAP_SUB);
        if (h2o_total <= 0.0)
            h2o_total = 0.0;
    }

    /*
 *  Determine total mass change due to evap/cond at soil
 */
    if (layer_count == 0)
        E_s_l = 0.0;
    else
    {
        if (layer_count == 2)
        {
            e_s_l = sati(T_s_l);
            T_bar = (T_g + T_s_l) / 2.0;
        }
        else
        {  /* layer_count == 1 */
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
        if (h2o_total > 0.0)
        {
            h2o_total += (E_s_l * VAP_SUB);
            if (h2o_total <= 0.0)
                h2o_total = 0.0;
        }
    }

    E_s = E_s_0 + E_s_l;

    /*      adj mass and depth for evap/cond        */

    if (layer_count > 0)
        _adj_snow(((E_s + (prev_h2o_tot - h2o_total)) / rho) / 2.0,
                  E_s);
}
/*
** NAME
**      _e_bal -- calculates point energy budget for 2-layer snowcover
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      int
**	_e_bal(void)
**
** DESCRIPTION
**      Calculates point energy budget for 2-layer snowcover.
**
** RETURN VALUE
**
**	TRUE	The calculations were completed.
**
**	FALSE	An error occured, and a message explaining the error has
**		been stored with the 'usrerr' routine.
**
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**
*/
int sno::_e_bal(void)
{
    if (snowcover)
    {

        /**	calculate energy xfer terms  **/

        /*      calculate net radiation */

        _net_rad();

        /*      calculate H & L_v_E  (and E as well)       */

        if (!_h_le())
            return 0;

        /*      calculate G & G_0(conduction/diffusion heat xfr)    */

        if (layer_count == 1)
        {
            G = g_soil(rho, T_s_0, T_g, z_s_0, z_g, P_a);
            G_0 = G;
        }
        else
        {  /*  layer_count == 2  */
            G = g_soil(rho, T_s_l, T_g, z_s_l, z_g, P_a);
            G_0 = g_snow(rho, rho, T_s_0, T_s_l, z_s_0, z_s_l,
                         P_a);
        }

        /*      calculate advection     */

        _advec();

        /*      sum E.B. terms  */

        /* surface energy budget */
        delta_Q_0 = R_n + H + L_v_E + G_0 + M;

        /* total snowpack energy budget */
        if (layer_count == 1)
            delta_Q = delta_Q_0;
        else  /* layer_count == 2 */
            delta_Q = delta_Q_0 + G - G_0;
    }
    else
    {
        R_n = 0.0;

        H = L_v_E = E = 0.0;

        G = G_0 = 0.0;

        M = 0.0;

        delta_Q = delta_Q_0 = 0.0;
    }

    return 1;
}


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


int sno::_do_tstep(
        TSTEP_REC *tstep)  /* timestep's record */
{
    /*
 *  A macro to update a time-weighted average for a quantity.
 *	avg		current average
 *	total_time	the time interval the current average applies to
 *	value		new value to be averaged in
 *	time_incr	the time interval the new value applies to
 */
#define TIME_AVG(avg, total_time, value, time_incr) \
        ( ((avg) * (total_time) + (value) * (time_incr)) \
        / ((total_time) + (time_incr)) )
    time_step = tstep->time_step;

    if (precip_now)
    {
        m_precip = precip_info[tstep->level].m_pp;
        m_rain = precip_info[tstep->level].m_rain;
        m_snow = precip_info[tstep->level].m_snow;
        z_snow = precip_info[tstep->level].z_snow;
    }

    h2o_total = 0.0;

    /*
     *  Is there a snowcover?
     */
    snowcover = (layer_count > 0);

    /*
     *  Calculate energy transfer terms
     */
    if (!_e_bal())
        return 0;

    /*
     *  Adjust mass and calculate runoff
     */
    _mass_bal();

    /*
     *  Update the averages for the energy terms and the totals for mass
     *  changes since the last output.
     */
    if (time_since_out > 0.0)
    {
        R_n_bar = TIME_AVG(R_n_bar, time_since_out,
                           R_n, time_step);
        H_bar = TIME_AVG(H_bar, time_since_out,
                         H, time_step);
        L_v_E_bar = TIME_AVG(L_v_E_bar, time_since_out,
                             L_v_E, time_step);
        G_bar = TIME_AVG(G_bar, time_since_out,
                         G, time_step);
        M_bar = TIME_AVG(M_bar, time_since_out,
                         M, time_step);
        delta_Q_bar = TIME_AVG(delta_Q_bar, time_since_out,
                               delta_Q, time_step);
        G_0_bar = TIME_AVG(G_0_bar, time_since_out,
                           G_0, time_step);
        delta_Q_0_bar = TIME_AVG(delta_Q_0_bar, time_since_out,
                                 delta_Q_0, time_step);

        E_s_sum += E_s;
        melt_sum += melt;
        ro_pred_sum += ro_predict;

        time_since_out += time_step;
    }
    else
    {
        R_n_bar = R_n;
        H_bar = H;
        L_v_E_bar = L_v_E;
        G_bar = G;
        M_bar = M;
        delta_Q_bar = delta_Q;
        G_0_bar = G_0;
        delta_Q_0_bar = delta_Q_0;

        E_s_sum = E_s;
        melt_sum = melt;
        ro_pred_sum = ro_predict;

        time_since_out = time_step;
    }

    /* increment time */
    current_time += time_step;

    if (tstep->output & WHOLE_TSTEP)
    {
        (*out_func)();
        if (!run_no_snow && (layer_count == 0))
            stop_no_snow = 1;
    }

    /*
     *  Update the model's input parameters
     */
    S_n += input_deltas[tstep->level].S_n;
    I_lw += input_deltas[tstep->level].I_lw;
    T_a += input_deltas[tstep->level].T_a;
    e_a += input_deltas[tstep->level].e_a;
    u += input_deltas[tstep->level].u;
    T_g += input_deltas[tstep->level].T_g;
    if (ro_data)
        ro += input_deltas[tstep->level].ro;

    return 1;
}

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
int sno::_divide_tstep(
        TSTEP_REC *tstep)    /* record of timestep to be divided */
{
    int next_level;    /* # of next level of timestep */
    TSTEP_REC *next_lvl_tstep;    /* info of next level of timestep */
    INPUT_REC *curr_lvl_deltas;    /* -> input-deltas of current level */
    INPUT_REC *next_lvl_deltas;    /* -> input-deltas of next level */
    PRECIP_REC *curr_lvl_precip;    /* -> precip data of current level */
    PRECIP_REC *next_lvl_precip;    /* -> precip data of next level */
    int i;            /* loop index */


    /*
     *  Fetch the record for the timestep at the next level.
     */
    next_level = tstep->level + 1;
    next_lvl_tstep = tstep_info + next_level;

    curr_lvl_deltas = input_deltas + tstep->level;
    next_lvl_deltas = input_deltas + next_level;

    curr_lvl_precip = precip_info + tstep->level;
    next_lvl_precip = precip_info + next_level;

    /*
     *  If this is the first time this new level has been used during
     *  the current data timestep, then calculate its input deltas
     *  and precipitation values.
     */
    if (!computed[next_level])
    {
        next_lvl_deltas->S_n = curr_lvl_deltas->S_n /
                               next_lvl_tstep->intervals;
        next_lvl_deltas->I_lw = curr_lvl_deltas->I_lw /
                                next_lvl_tstep->intervals;
        next_lvl_deltas->T_a = curr_lvl_deltas->T_a /
                               next_lvl_tstep->intervals;
        next_lvl_deltas->e_a = curr_lvl_deltas->e_a /
                               next_lvl_tstep->intervals;
        next_lvl_deltas->u = curr_lvl_deltas->u /
                             next_lvl_tstep->intervals;
        next_lvl_deltas->T_g = curr_lvl_deltas->T_g /
                               next_lvl_tstep->intervals;
        if (ro_data)
            next_lvl_deltas->ro = curr_lvl_deltas->ro /
                                  next_lvl_tstep->intervals;

        if (precip_now)
        {
            next_lvl_precip->m_pp = curr_lvl_precip->m_pp /
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

//    if(!_do_tstep(next_lvl_tstep))
//        return 0;
//    else
//        return 1;

    /*
     *  For each the new smaller timestep, either subdivide them if
     *  below their mass threshold, or run the model for them.
     */
    for (i = 0; (i < next_lvl_tstep->intervals) ; i++) //&& !stop_no_snow  <-- not used because we aren't using output funcs
    {
        if ((next_level != SMALL_TSTEP) &&
            _below_thold(next_lvl_tstep->threshold))
        {
            if (!_divide_tstep(next_lvl_tstep))
                return 0;
        }
        else
        {
            if (!_do_tstep(next_lvl_tstep))
                return 0;
        }
    }


    /*
     *  Output if this timestep is divided?
     */
    if (tstep->output & DIVIDED_TSTEP)
    {
        (*out_func)();
        if (!run_no_snow && (layer_count == 0))
            stop_no_snow = 1;
    }

    return 1;
}

/*
** NAME
**      _cold_content -- calculates cold content for a layer
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      double
**	_cold_content(
**	    double  temp,		|* temperature of layer *|
**	    double  mass)		|* specific mass of layer *|
**
** DESCRIPTION
**      This routine calculates the cold content for a layer (i.e., the
**	energy required to bring its temperature to freezing) from the
**	layer's temperature and specific mass.
**
** RETURN VALUE
**	The layer's cold content.
*/

double sno::_cold_content(
        double temp,        /* temperature of layer */
        double mass)        /* specific mass of layer */
{
    if (temp < FREEZE)
        return heat_stor(CP_ICE(temp), mass, (temp - FREEZE));
    else
        return 0.0;
}

/*
** NAME
**      _calc_layers -- determine # of layers in snowcover and their depths
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_calc_layers(void);
**
** DESCRIPTION
**      This routine determines the # of layers in the snowcover based its
**	depth and mass.  Usually, there are are 2 layers: the surface (active)
**	and the lower layer.  The depth of the surface layer is set to the
**	maximum depth for the surface layer (variable "max_z_s_0").  The
**	remaining depth constitutes the lower layer.  The routine checks
**	to see if the mass of this lower layer is above the minium threshold
**	(i.e., the mass threshold for the small run timestep).  If not,
**	the surface layer is the whole snowcover, and there's no lower
**	layer.
**
** GLOBAL VARIABLES READ
**	m_s
**	max_z_s_0
**	rho
**	tstep_info
**	z_s
**
** GLOBAL VARIABLES MODIFIED
**	layer_count
**	z_s
**	z_s_0
**	z_s_l
*/
void sno::_calc_layers(void)
{
    if (m_s <= tstep_info[SMALL_TSTEP].threshold)
    {
        /*
         *  Less than minimum layer mass, so treat as no snowcover.
         */
        layer_count = 0;
        z_s = z_s_0 = z_s_l = 0.0;
    }
    else if (z_s < max_z_s_0)
    {
        /*
         *  Not enough depth for surface layer and the lower layer,
         *  so just 1 layer: surface layer.
         */
        layer_count = 1;
        z_s_0 = z_s;
        z_s_l = 0.0;
    }
    else
    {
        /*
         *  Enough depth for both layers.
         */
        layer_count = 2;
        z_s_0 = max_z_s_0;
        z_s_l = z_s - z_s_0;

        /*
         *  However, make sure there's enough MASS for the lower
         *  layer.  If not, then there's only 1 layer.
         */
        if (z_s_l * rho < tstep_info[SMALL_TSTEP].threshold)
        {
            layer_count = 1;
            z_s_0 = z_s;
            z_s_l = 0.0;
        }
    }
}

/*
** NAME
**      _below_thold -- is a layer's mass below a threshold ?
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      int
**	_below_thold(
**	    double  threshold)	|* current timestep's threshold for a
**				   layer's mass *|
**
** DESCRIPTION
**      This routine determines if any individual layer's mass is below
**	a given threshold for the current timestep.
**
** RETURN VALUE
**      	1	A layer's mass is less than the threshold.
**
**		0	All layers' masses are greater than the threshold.
**
** GLOBAL VARIABLES READ
**	layer_count
**	m_s
**	m_s_0
**	m_s_l
**
** GLOBAL VARIABLES MODIFIED
*/
int sno::_below_thold(
        double threshold)    /* current timestep's threshold for a
				   layer's mass */
{
    if (layer_count == 0)
        return 0;
    if (layer_count == 1)
        return (m_s < threshold);
    else  /* layer_count == 2 */
        return ((m_s_0 < threshold) || (m_s_l < threshold));
}

/*
** NAME
        **      _advec -- calculates advected energy at a point
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
        **	_advec(void)
**
** DESCRIPTION
        **      This routine calculates the advected energy for a 2-layer snowcover
**	if there's precipitation for the current timestep.
**
** GLOBAL VARIABLES READ
        **	m_rain
**	m_snow
        **	precip_now
**	T_rain
        **	T_s_0
**	T_snow
        **	time_step
**
** GLOBAL VARIABLES MODIFIED
        **	M
*/
void sno::_advec(void)
{
    if (precip_now)
    {
        M = (heat_stor(CP_WATER(T_rain), m_rain, (T_rain - T_s_0)) +
             heat_stor(CP_ICE(T_snow), m_snow, (T_snow - T_s_0)))
            / time_step;
    }
    else
    {
        M = 0.0;
    }

    if (M < -1000)
    {
        M=0;
        //LOG_DEBUG << "WUUUUUUUT";
    }
}

/** NAME
**      _adj_snow -- adjust the snowcover with changes in depth and/or mass
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_adj_snow(
**	    double delta_z_s,	|* change in snowcover's depth *|
**	    double delta_m_s)	|* change is snowcover's mass *|
**
** DESCRIPTION
**      This routine adjusts the snowcover for a change in its depth or
**	its mass or both.  The snowcover's density is updated.  If there
**	is a change in the snowcover's depth, the # of layers is recomputed.
**	If there's just a change in the snowcover's mass with no change in
**	its depth, then just the specific masses for the layers are updated.
**
**	The routine ensures that the snowcover's density does NOT exceed
**	a maximum density (currently 750 kg/m^3).  If the adjustments to
**	the snowcover, for some reason, lead to an excessive density, the
**	density is clipped at the maximum, and the depth re-adjusted
**	accordingly.
**
** GLOBAL VARIABLES READ
**
** GLOBAL VARIABLES MODIFIED
**	m_s
**	rho
**	z_s
**
**	(and those variables modified by "_adj_layers" and "_layer_mass")
*/
void sno::_adj_snow(
        double delta_z_s,    /* change in snowcover's depth */
        double delta_m_s)    /* change is snowcover's mass */
{

    /*
     *  Update depth, mass, and then recompute density.
     */
    z_s += delta_z_s;
    m_s += delta_m_s;


    if(m_s < 0 || z_s < 0)
    {
        m_s = cc_s = 0.0;
        m_s_0 = cc_s_0 = 0.0;
    }
    if (z_s != 0.0)
    {
        rho = m_s / z_s;
    } else
    {
        rho = 0;
    }

    /*
     *  Clip density at maxium density if necessary.
     */
    if (rho > MAX_SNOW_DENSITY)
    {
        rho = MAX_SNOW_DENSITY;
        z_s = m_s / rho;
        _adj_layers();
    }
    else
    {
        /*
         *  If a change in depth, adjust the layers' depths and masses.
         */
        if (delta_z_s != 0.0)
            _adj_layers();
        else
            /*
              *  Just a change in the snowcover's mass, so update the
              *  layers' masses.
              */
            _layer_mass();
    }
}

/*
** NAME
**      _adj_layers -- adjust the layers because of new snowcover depth
**
** SYNOPSIS
**      #include "_snobal.h"
**
**      void
**	_adj_layers(void);
**
** DESCRIPTION
**      This routine adjusts the layers of the snowcover because the
**	snowcover's depth has changed.  It is assumed that the snowcover's
**	density has already been updated.  The # of layers are recomputed
**	based on the overall snowcover depth.  Their depths and masses
**	are updated as well.  If a layer has been created due to an
**	increase in the snowcover's depth, its temperature and cold content
**	are initialized.
**
** GLOBAL VARIABLES READ
**	layer_count
**
** GLOBAL VARIABLES MODIFIED
**	cc_s
**	cc_s_0
**	cc_s_l
**	h2o
**	h2o_max
**	h2o_total
**	h2o_vol
**	m_s
**	m_s_0
**	m_s_l
**	rho
**	T_s
**	T_s_0
**	T_s_l
**
**	(and those variables modified by "_calc_layer" and "_layer_mass")
*/
void sno::_adj_layers(void)
{
    int prev_layer_count;    /* previous # of layers, if change in depth */

    /*
     *  Recompute then number of layers and see if there's been
     *  a change in the # of layers.  Note:  since this routine
     *  is called to adjust an existing snowcover, the current # of
     *  layers must be either 1 or 2 while the new # of layers may
     *  either be 0, 1 or 2.
     *
     *	current #	new #
     *	of layers	of layers
     *
     *	   1	   -->	   0
     *	   1	   -->	   1	(no change)
     *	   1	   -->	   2
     *	   2	   -->	   0
     *	   2	   -->	   1
     *	   2	   -->	   2	(no change)
     */
    prev_layer_count = layer_count;  /* must be > 0 */
    _calc_layers();

    if (layer_count == 0)
    {
        /*
         *  1 or 2 layers --> 0 layers
         */
        rho = 0.0;

        /*
         *  If mass > 0, then it must be below threshold.
         *  So turn this little bit of mass into water.
         */
        if (m_s > 0.0)
            h2o_total += m_s;

        m_s = cc_s = 0.0;
        m_s_0 = cc_s_0 = 0.0;

        /*
         *  Note: Snow temperatures are set to MIN_SNOW_TEMP
         *	  (as degrees K) instead of 0 K to keep quantization
         *	  range in output image smaller.
         */
        T_s = T_s_0 = MIN_SNOW_TEMP + FREEZE;

        if (prev_layer_count == 2)
        {
            m_s_l = cc_s_l = 0.0;
            T_s_l = MIN_SNOW_TEMP + FREEZE;
        }
        h2o_vol = h2o = h2o_max = h2o_sat = 0.0;
    }

    else
    {
        _layer_mass();

        if ((prev_layer_count == 1) && (layer_count == 2))
        {
            /*
             *  1 layer --> 2 layers, add lower layer
             */
            T_s_l = T_s;
            cc_s_l = _cold_content(T_s_l, m_s_l);
        }

        else if ((prev_layer_count == 2) && (layer_count == 1))
        {
            /*
             *  2 layers --> 1 layer, remove lower layer
             */
            T_s_l = MIN_SNOW_TEMP + FREEZE;
            cc_s_l = 0.0;
        }
    }
}
