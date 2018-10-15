/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "snomacros.h"
#include <cmath>
#include "logger.hpp"
#include "exception.hpp"

typedef struct
{
    double m_pp;
    /* total precipitation mass (kg/m^2) */
    double m_rain;
    /* mass of rain in precip (kg/m^2) */
    double m_snow;
    /*  "   "  snow "     "   (kg/m^2) */
    double z_snow;    /* depth of snow in   "   (m) */
} PRECIP_REC;
typedef struct
{
    double S_n;
    /* net solar radiation (W/m^2) */
    double I_lw;
    /* incoming longwave (thermal) rad (W/m^2) */
    double T_a;
    /* air temp (C) */
    double e_a;
    /* vapor pressure (Pa) */
    double u;
    /* wind speed (m/sec) */
    double T_g;
    /* soil temp at depth z_g (C) */
    double ro;    /* measured runoff (m/sec) */
} INPUT_REC;
typedef struct
{
    int level;    /* timestep's level */
#define    DATA_TSTEP       0
#define    NORMAL_TSTEP       1
#define    MEDIUM_TSTEP       2
#define    SMALL_TSTEP       3

    double time_step;
    /* length of timestep (seconds) */
    int intervals;
    /* # of these timestep that are in
                       the previous-level's timestep
                       (not used for level 0: data tstep) */
    double threshold;
    /* mass threshold for a layer to use
                       this timestep
                       (not used for level 0: data tstep) */
    int output;    /* flags whether or not to call output
					   function for timestep */
#define WHOLE_TSTEP      0x1        /* output when tstep is not divided */
#define DIVIDED_TSTEP      0x2        /* output when timestep is divided */

} TSTEP_REC;

class sno
{
public:
    /*
     * Functions
     */

    double ssxfr(
            double k1,    /* layer 1's thermal conductivity (J / (m K sec))  */
            double k2,    /* layer 2's    "         "                        */
            double t1,    /* layer 1's average layer temperature (K)	   */
            double t2,    /* layer 2's    "      "        "         	   */
            double d1,     /* layer 1's thickness (m)			   */
            double d2);     /* layer 2's    "       "			   */

    double satw(double tk);        /* air temperature (K)		*/
    double sati(double tk);        /* air temperature (K)	*/
    double new_tsno(
            double spm,    /* layer's specific mass (kg/m^2) 	 */
            double t0,    /* layer's last temperature (K) 	 */
            double ccon);    /* layer's adjusted cold content (J/m^2) */
    void init_snow(void);
    int hle1(
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
            double *e);    /* mass flux (+ to surf) (kg/m^2/s)	*/

    double heat_stor(
            double cp,    /* specific heat of layer (J/kg K) */
            double spm,    /* layer specific mass (kg/m^2)    */
            double tdif);    /* temperature change (K)          */

    double g_soil(
            double rho,    /* snow layer's density (kg/m^3)	     */
            double tsno,    /* snow layer's temperature (K)		     */
            double tg,    /* soil temperature (K)			     */
            double ds,    /* snow layer's thickness (m)		     */
            double dg,    /* dpeth of soil temperature measurement (m) */
            double pa);    /* air pressure (Pa)			     */
    double g_snow(
            double rho1,    /* upper snow layer's density (kg/m^3)	*/
            double rho2,    /* lower  "     "        "    (kg/m^3)	*/
            double ts1,    /* upper snow layer's temperature (K)	*/
            double ts2,    /* lower  "     "         "       (K)	*/
            double ds1,    /* upper snow layer's thickness (m)	*/
            double ds2,    /* lower  "     "         "     (m)	*/
            double pa);    /* air pressure (Pa)      			*/
    double efcon(
            double k,    /* layer thermal conductivity (J/(m K sec)) */
            double t,    /* layer temperature (K)		    */
            double p);    /* air pressure (Pa)  			    */
    int do_data_tstep(void);

    void _time_compact_ori(void);
    void _time_compact(void);
    void    _snowmelt(void);

    void _runoff(void);

    void _new_density(void);
    void _precip(void);
    void _net_rad(void);
    void _mass_bal(void);
    void _layer_mass(void);
    int _h_le(void);
    void _h2o_compact(void);
    void _evap_cond(void);
    int _e_bal(void);
    int _do_tstep(
            TSTEP_REC *tstep);  /* timestep's record */
    double _cold_content(
            double temp,        /* temperature of layer */
            double mass);        /* specific mass of layer */
    int _divide_tstep(
            TSTEP_REC *tstep);    /* record of timestep to be divided */
    void _calc_layers(void);
    int _below_thold(
            double threshold);    /* current timestep's threshold for a
				   layer's mass */
    void _advec(void);
    void _adj_snow(
            double delta_z_s,    /* change in snowcover's depth */
            double delta_m_s);    /* change is snowcover's mass */
    void _adj_layers(void);
    double psi(
            double zeta,        /* z/lo				*/
            int code);        /* which psi function? (see above) */


    /**
     * Debug
     *
     */
    size_t _debug_id;
/*
 *  Global variables that are used to communicate with the snobal library
 *  routines.  These are public, i.e., these can be accessed from outside
 *  the library.
 */

    INPUT_REC input_deltas[4];
    /* deltas for climate-input parameters
                       over each timestep */

    PRECIP_REC precip_info[4];
    /* array of precip info adjusted for
                       each timestep */
    int errno;

    //thermal conductivity of wet sand (J/(m sec K))
    double KT_WETSAND;

    int computed[4];
    /* array of flags for each timestep;
                   TRUE if computed values for input
                   deltas and precip arrays */

    int isothermal;
    /* melting? */
    int snowcover;      /* snow on gnd at start of current timestep? */

/*   variables that control model execution   */

    int run_no_snow;
    /* continue model even if snow disappears? */
    int stop_no_snow;

    /* stopped model because no snow left? */

    void    (*out_func)(void);    /* -> output function */


/*   constant model parameters  */

    double max_z_s_0;
    /* maximum active layer thickness (m) */
    double max_h2o_vol;    /* max liquid h2o content as volume ratio:
				     V_water/(V_snow - V_ice) (unitless) */


/*   time step information */

    TSTEP_REC tstep_info[4];
    /* array of info for each timestep:
                          0 : data timestep
                          1 : normal run timestep
                          2 : medium  "     "
                          3 : small   "     "
                    */

    double time_step;
    /* length current timestep (sec) */
    double current_time;
    /* start time of current time step (sec) */
    double time_since_out;    /* time since last output record (sec) */


/*   snowpack information   */

    int layer_count;
    /* number of layers in snowcover: 0, 1, or 2 */
    double z_s;
    /* total snowcover thickness (m) */
    double z_s_0;
    /* active layer depth (m) */
    double z_s_l;
    /* lower layer depth (m) */
    double rho;
    /* average snowcover density (kg/m^3) */
    double m_s;
    /* snowcover's specific mass (kg/m^2) */
    double m_s_0;
    /* active layer specific mass (kg/m^2) */
    double m_s_l;
    /* lower layer specific mass (kg/m^2) */
    double T_s;
    /* average snowcover temp (K) */
    double T_s_0;
    /* active snow layer temp (K) */
    double T_s_l;
    /* lower layer temp (C) */
    double cc_s;
    /* snowcover's cold content (J/m^2) */
    double cc_s_0;
    /* active layer cold content (J/m^2) */
    double cc_s_l;
    /* lower layer cold content (J/m^2) */
    double h2o_sat;
    /* % of liquid H2O saturation (relative water
                 content, i.e., ratio of water in snowcover
                 to water that snowcover could hold at
                 saturation) */
    double h2o_vol;
    /* liquid h2o content as volume ratio:
                 V_water/(V_snow - V_ice) (unitless) */
    double h2o;
    /* liquid h2o content as specific mass
             (kg/m^2) */
    double h2o_max;
    /* max liquid h2o content as specific mass
                 (kg/m^2) */
    double h2o_total;      /* total liquid h2o: includes h2o in snowcover,
				     melt, and rainfall (kg/m^2) */


/*   climate-data input records   */

    int ro_data;
    /* runoff data? */

    INPUT_REC input_rec1;
    /* input data for start of data timestep */
    INPUT_REC input_rec2;
    /*   "     "   "  end   "   "      "     */

/*   climate-data input values for the current run timestep */

    double S_n;
    /* net solar radiation (W/m^2) */
    double I_lw;
    /* incoming longwave (thermal) rad (W/m^2) */
    double T_a;
    /* air temp (C) */
    double e_a;
    /* vapor pressure (Pa) */
    double u;
    /* wind speed (m/sec) */
    double T_g;
    /* soil temp at depth z_g (C) */
    double ro;             /* measured runoff (m/sec) */


/*   other climate input   */

    double P_a;            /* air pressure (Pa) */


/*   measurement heights/depths   */

    int relative_hts;
    /* TRUE if measurements heights, z_T
                   and z_u, are relative to snow
                   surface; FALSE if they are
                   absolute heights above the ground */
    double z_g;
    /* depth of soil temp meas (m) */
    double z_u;
    /* height of wind measurement (m) */
    double z_T;
    /* height of air temp & vapor pressure
           measurement (m) */
    double z_0;            /* roughness length */

/*   precipitation info for the current DATA timestep    */

    int precip_now;
    /* precipitation occur for current timestep? */
    double m_pp;
    /* specific mass of total precip (kg/m^2) */
    double percent_snow;
    /* % of total mass that's snow (0 to 1.0) */
    double rho_snow;
    /* density of snowfall (kg/m^3) */
    double T_pp;
    /* precip temp (C) */
    double T_rain;
    /* rain's temp (K) */
    double T_snow;
    /* snowfall's temp (K) */
    double h2o_sat_snow;   /* snowfall's % of liquid H2O saturation */

/*   precipitation info adjusted for current run timestep   */

    double m_precip;
    /* specific mass of total precip (kg/m^2) */
    double m_rain;
    /*    "      "   of rain in precip (kg/m^2) */
    double m_snow;
    /*    "      "   "  snow "    "    (kg/m^2) */
    double z_snow;        /* depth of snow in precip (m) */


/*   energy balance info for current timestep        */

    double R_n;
    /* net allwave radiation (W/m^2) */
    double H;
    /* sensible heat xfr (W/m^2) */
    double L_v_E;
    /* latent heat xfr (W/m^2) */
    double G;
    /* heat xfr by conduction & diffusion from soil
           to snowcover (W/m^2) */
    double G_0;
    /* heat xfr by conduction & diffusion from soil
             or lower layer to active layer (W/m^2) */
    double M;
    /* advected heat from precip (W/m^2) */
    double delta_Q;
    /* change in snowcover's energy (W/m^2) */
    double delta_Q_0;      /* change in active layer's energy (W/m^2) */

/*   averages of energy balance vars since last output record   */

    double R_n_bar;
    double H_bar;
    double L_v_E_bar;
    double G_bar;
    double G_0_bar;
    double M_bar;
    double delta_Q_bar;
    double delta_Q_0_bar;


/*   mass balance vars for current timestep        */

    double melt;
    /* specific melt (kg/m^2 or m) */
    double E;
    /* mass flux by evap into air from active
                 layer (kg/m^2/s) */
    double E_s;
    /* mass of evap into air & soil from snowcover
                 (kg/m^2) */
    double ro_predict;     /* predicted specific runoff (m/sec) */

/*   sums of mass balance vars since last output record   */

    double melt_sum;
    double E_s_sum;
    double ro_pred_sum;

};

