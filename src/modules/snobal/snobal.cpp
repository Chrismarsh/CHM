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

#include "snobal.hpp"
REGISTER_MODULE_CPP(snobal);

snobal::snobal(config_file cfg)
        : module_base("snobal", parallel::data, cfg)
{
    depends("frac_precip_snow");
    depends("iswr");
    depends("rh");
    depends("t");
    depends("U_2m_above_srf");
    depends("p");
    depends("ilwr");

    // Optional subcanopy variables if a canopy module is included (used if exist)
    optional("frac_precip_snow_subcanopy");
    optional("iswr_subcanopy");
    optional("rh_subcanopy");
    optional("ta_subcanopy");
    optional("p_subcanopy");
    optional("ilwr_subcanopy");

    optional("drift_depth");
    optional("drift_mass");

    depends("snow_albedo");

    optional("T_g");

    // Optional avalanche variables
    optional("delta_avalanche_snowdepth");
    optional("delta_avalanche_mass");

    provides("swe");
    provides("snowmelt_int");
    provides("R_n");
    provides("H");
    provides("E");
    provides("G");
    provides("M");
    provides("dQ");
    provides("cc");
    provides("T_s");
    provides("T_s_0");
    provides("T_s_l");
    provides("dead");
    provides("iswr_net");
    provides("isothermal");
    provides("ilwr_out");
    provides("sum_snowpack_runoff");
    provides("sum_melt");
    provides("snowdepthavg");
    provides("snowdepthavg_vert");


}

void snobal::init(mesh& domain)
{
    ompException oe;

    drift_density = cfg.get("drift_density",300.);
    const_T_g = cfg.get("const_T_g",-4.0);

    //store all of snobals global_param variables from this timestep to be used as ICs for the next timestep
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
      oe.Run([&]
	     {
	       auto face = domain->face(i);

	       snodata* g = face->make_module_data<snodata>(ID);
	       g->sum_runoff = 0;
	       g->sum_melt = 0;
	       auto* sbal = &(g->data);
	       g->dead=0;
	       g->delta_avalanche_snowdepth=0;
	       g->delta_avalanche_swe=0;
	       /**
		* Snopack config
		*/
	       sbal->param_snow_compaction = cfg.get("param_snow_compaction",1); // new param is the default

	       sbal->h2o_sat = .3;
	       sbal->layer_count = 0;
	       sbal->m_s = 0.;
	       sbal->m_s_0 = 0.;
	       sbal->m_s_l = 0.;
	       sbal->max_h2o_vol = cfg.get("max_h2o_vol",.0001);//0.0001
	       sbal->rho = 0.;
	       sbal->T_s = -75. + FREEZE;
	       sbal->T_s_0 = -75. + FREEZE; //assuming no snow
	       sbal->T_s_l = -75. + FREEZE;
	       sbal->z_s = 0.;

	       sbal->KT_WETSAND = cfg.get("kt_wetsand",0.08);

	       sbal->ro_data = 0;

	       sbal->max_z_s_0 = cfg.get("max_active_layer",.1);
	       sbal->h2o_total = 0;
	       sbal->isothermal = 0;

	       /// Heights
	       sbal->z_0 = cfg.get("z_0",0.001);
	       sbal->z_T = cfg.get("z_T",2.6);
	       sbal->z_u = cfg.get("z_u",2.96);
	       sbal->z_g = cfg.get("z_g",0.1);
	       sbal->relative_hts = 1; // True (1) -- relative to the snow surface via scale_wind_speed which takes into account snowdepth.

	       sbal->R_n_bar = 0.0;
	       sbal->H_bar = 0.0;
	       sbal->L_v_E_bar = 0.0;
	       sbal->G_bar = 0.0;
	       sbal->M_bar = 0.0;
	       sbal->delta_Q_bar = 0.0;
	       sbal->E_s_sum = 0.0;
	       sbal->melt_sum = 0.0;
	       sbal->ro_pred_sum = 0.0;

	       sbal->time_since_out = 0;
	       sbal->current_time = 0;
	       sbal->run_no_snow = 1;
	       sbal->stop_no_snow = 1;
	       sbal->snowcover = 0;
	       sbal->precip_now = 0;

	       //init the step_info struct
	       sbal->tstep_info[DATA_TSTEP].level = DATA_TSTEP;
	       sbal->tstep_info[DATA_TSTEP].time_step = global_param->dt();
	       sbal->tstep_info[DATA_TSTEP].intervals = 0;
	       sbal->tstep_info[DATA_TSTEP].threshold = 20;
	       sbal->tstep_info[DATA_TSTEP].output = 0;


	       sbal->tstep_info[NORMAL_TSTEP].level = NORMAL_TSTEP;
	       sbal->tstep_info[NORMAL_TSTEP].time_step = global_param->dt();
	       sbal->tstep_info[NORMAL_TSTEP].intervals = sbal->tstep_info[DATA_TSTEP].time_step /
		 sbal->tstep_info[NORMAL_TSTEP].time_step;;
	       sbal->tstep_info[NORMAL_TSTEP].threshold = 20;
	       sbal->tstep_info[NORMAL_TSTEP].output = 0;


	       sbal->tstep_info[MEDIUM_TSTEP].level = MEDIUM_TSTEP;
	       sbal->tstep_info[MEDIUM_TSTEP].time_step = global_param->dt()/4;
	       sbal->tstep_info[MEDIUM_TSTEP].intervals = sbal->tstep_info[NORMAL_TSTEP].time_step /
		 sbal->tstep_info[MEDIUM_TSTEP].time_step;
	       sbal->tstep_info[MEDIUM_TSTEP].threshold = 10;
	       sbal->tstep_info[MEDIUM_TSTEP].output = 0;


	       sbal->tstep_info[SMALL_TSTEP].level = SMALL_TSTEP;
	       sbal->tstep_info[SMALL_TSTEP].time_step = global_param->dt()/ 100;// 60;
	       sbal->tstep_info[SMALL_TSTEP].intervals = sbal->tstep_info[MEDIUM_TSTEP].time_step /
		 sbal->tstep_info[SMALL_TSTEP].time_step;
	       sbal->tstep_info[SMALL_TSTEP].threshold = 0.2;
	       sbal->tstep_info[SMALL_TSTEP].output = 0;


	       ////////
	       if (face->has_initial_condition("swe"))
	       {
		   if( !is_nan(face->get_initial_condition("swe")))
		   {
		       sbal->rho = cfg.get("IC_rho",300.);
		       sbal->T_s =  -10 + FREEZE;
		       sbal->T_s_0 = -15. + FREEZE; //assuming no snow
		       sbal->T_s_l = -15. + FREEZE;
		       sbal->z_s = face->get_initial_condition("swe") / sbal->rho;

		   }
	       }

	       sbal->init_snow();

	       //in point mode, the entire mesh still exists, but no timeseries has been allocated for the faces
	       //thus this segfaults. This should be fixed

	       if (face->has_initial_condition("swe") &&  !is_nan(face->get_initial_condition("swe")))
	       {
		   (*face)["swe"_s]= sbal->m_s;
		   (*face)["R_n"_s]= sbal->R_n;
		   (*face)["H"_s]= sbal->H;
		   (*face)["E"_s]= sbal->L_v_E;
		   (*face)["G"_s]= sbal->G;
		   (*face)["M"_s]= sbal->M;
		   (*face)["dQ"_s]= sbal->delta_Q;
		   (*face)["cc"_s]= sbal->cc_s;
		   (*face)["T_s"_s]= sbal->T_s;
		   (*face)["T_s_0"_s]= sbal->T_s_0;
		   (*face)["T_s_l"_s]= sbal->T_s_l;
		   (*face)["iswr_net"_s]= sbal->S_n;
		   (*face)["isothermal"_s]= sbal->isothermal;
		   (*face)["ilwr_out"_s]= sbal->R_n - sbal->S_n - sbal->I_lw;
		   (*face)["snowmelt_int"_s]= 0;
		   (*face)["sum_melt"_s]= g->sum_melt;
		   (*face)["sum_snowpack_runoff"_s]= g->sum_runoff;

		   (*face)["snowdepthavg"_s]= sbal->z_s;
	       }
	     });
    }
    oe.Rethrow();
}
snobal::~snobal()
{

}

void snobal::run(mesh_elem &face)
{
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }

    //debugging
    auto id = face->cell_id;

    bool run=false;


    auto hour = global_param->hour();
    auto day = global_param->day();
    auto month = global_param->month();


    //get the previous timesteps data out of the global_param store.
    snodata* g = face->get_module_data<snodata>(ID);
    auto* sbal = &(g->data);

    sbal->_debug_id = id;

    sbal->P_a = mio::Atmosphere::stdAirPressure( face->get_z());

    double albedo = (*face)["snow_albedo"_s];

    // Optional inputs if there is a canopy or not
    double ilwr;
    if(has_optional("ilwr_subcanopy")) {
        ilwr = (*face)["ilwr_subcanopy"_s];
    } else {
        ilwr = (*face)["ilwr"_s];
    }

    // Optional inputs if there is a canopy or not
    double rh;
    if(has_optional("rh_subcanopy")) {
        rh = (*face)["rh_subcanopy"_s];
    } else {
        rh = (*face)["rh"_s];
    }

    // Optional inputs if there is a canopy or not
    double t;
    if(has_optional("ta_subcanopy")) {
        t = (*face)["ta_subcanopy"_s];
    } else {
        t = (*face)["t"_s];
    }

    double ea = mio::Atmosphere::vaporSaturationPressure(t+273.15)  * rh/100.;

    // Optional inputs if there is a canopy or not
    if(has_optional("iswr_subcanopy")) {
        sbal->input_rec2.S_n = (1.0-albedo) * (*face)["iswr_subcanopy"_s];
    } else {
        sbal->input_rec2.S_n = (1.0-albedo) * (*face)["iswr"_s];
    }
    sbal->input_rec2.I_lw = ilwr;
    sbal->input_rec2.T_a = t+FREEZE;
    sbal->input_rec2.e_a = ea;

    // Optional inputs if there is a canopy or not
    sbal->input_rec2.u = (*face)["U_2m_above_srf"_s];
    sbal->input_rec2.u = std::max(sbal->input_rec2.u,1.0);

    if(has_optional("T_g"))
        sbal->input_rec2.T_g = (*face)["T_g"_s] + 273.15;
    else
        sbal->input_rec2.T_g = const_T_g+FREEZE;


    sbal->input_rec2.ro = 0.;

    if(global_param->first_time_step)
    {
        sbal->input_rec1.S_n = sbal->input_rec2.S_n;
        sbal->input_rec1.I_lw =  sbal->input_rec2.I_lw;
        sbal->input_rec1.T_a =  sbal->input_rec2.T_a;
        sbal->input_rec1.e_a =  sbal->input_rec2.e_a;
        sbal->input_rec1.u = sbal->input_rec2.u;
        sbal->input_rec1.T_g =  sbal->input_rec2.T_g; //TODO: FIx this with a gflux estimate
        sbal->input_rec1.ro =  sbal->input_rec2.ro;
    }

    // Optional inputs if there is a canopy or not
    double p;
    if(has_optional("p_subcanopy")) {
        p = (*face)["p_subcanopy"_s];
    } else {
        p = (*face)["p"_s];
    }

    if(p >= 0.00025) //0.25mm swe
    {
        sbal->precip_now = 1;
        sbal->m_pp = p;

        // Optional inputs if there is a canopy or not
        if(has_optional("frac_precip_snow_subcanopy")) {
            sbal->percent_snow = (*face)["frac_precip_snow_subcanopy"_s];
        } else {
            sbal->percent_snow = (*face)["frac_precip_snow"_s];
        }
        sbal->rho_snow = 100.; //http://ccc.atmos.colostate.edu/pdfs/SnowDensity_BAMS.pdf
        sbal->T_pp = t; //actually in C unlike everything else in the model!!  //+FREEZE;//std::min(t+FREEZE,0.0);
        sbal->stop_no_snow=0;
    }
    else
    {
        sbal->precip_now = 0;
        sbal->stop_no_snow=0;
    }

    if(has_optional("drift_mass"))
    {

        double mass = (*face)["drift_mass"_s];
        mass = is_nan(mass) ? 0 : mass;
        //m_s is kg/m^2 and mass is kg/m^2
        //negative = mass removal
//        if(mass < 0 && (sbal->m_s+mass ) < 0 ) // are we about to remove more mass than what exists???
//            mass = -sbal->m_s; //cap it to remove no more than available mass

        sbal->_adj_snow(mass / drift_density, mass);
    }

    // If snow avalanche variables are available
    bool snow_slide = false;
    if(has_optional("delta_avalanche_snowdepth")) {
        g->delta_avalanche_snowdepth = (*face)["delta_avalanche_snowdepth"_s];
    }
    if(has_optional("delta_avalanche_mass")) {
        g->delta_avalanche_swe = (*face)["delta_avalanche_mass"_s];
        snow_slide = true;
    }

    // Redistribute snow (if snow_slide is used)
    if(snow_slide) {
        // _adj_snow(depth change (m), swe change (kg/m^2))
        // Convert change in volume and mass back to depth and mass per area, respectivly.
        // Assumes snow depth is uniform across triangle
        double area = face->get_area(); // area of current triangle (m^2)
        double d_depth = g->delta_avalanche_snowdepth / area; // m^3 / m^2 = m
        double d_mass  = g->delta_avalanche_swe / area * 1000; // m^3 / m^2 * 1000 kg/m^3 = kg/m^2
        sbal->_adj_snow(d_depth,d_mass);
    }


    if(g->dead == 1)
    {
        sbal->init_snow();
        g->dead = 0;
    }

    double swe1 = sbal->m_s;
    try
    {
        sbal->do_data_tstep();
    }catch(module_error& e)
    {
        g->dead=1;
        LOG_DEBUG << boost::diagnostic_information(e);
        auto details = "("+std::to_string(face->center().x()) + "," + std::to_string(face->center().y())+","+std::to_string(face->center().z())+") ID = " + std::to_string(face->cell_id);
//        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Snobal died. Triangle center = "+details));
    }



    double swe_diff = swe1 - sbal->m_s;
    swe_diff = swe_diff > 0. ? swe_diff : 0;
    g->sum_runoff += sbal->ro_predict;
    g->sum_melt += swe_diff;

    double sd_ver = sbal->z_s/std::max(0.001,cos(face->slope()));

    (*face)["dead"_s]=g->dead;

    (*face)["swe"_s]=sbal->m_s;

    (*face)["R_n"_s]=sbal->R_n;
    (*face)["H"_s]=sbal->H;
    (*face)["E"_s]=sbal->L_v_E;
    (*face)["G"_s]=sbal->G;
    (*face)["M"_s]=sbal->M;
    (*face)["dQ"_s]=sbal->delta_Q;
    (*face)["cc"_s]=sbal->cc_s;
    (*face)["T_s"_s]=sbal->T_s;
    (*face)["T_s_0"_s]=sbal->T_s_0;
    (*face)["T_s_l"_s]=sbal->T_s_l;
    (*face)["iswr_net"_s]=sbal->S_n;
    (*face)["isothermal"_s]=sbal->isothermal;
    (*face)["ilwr_out"_s]= sbal->R_n - sbal->S_n - sbal->I_lw;
//    (*face)["snowmelt_int"_s]=sbal->ro_predict;

    (*face)["snowmelt_int"_s]=swe_diff;
    (*face)["sum_melt"_s]=g->sum_melt;
    (*face)["sum_snowpack_runoff"_s]=g->sum_runoff;

    (*face)["snowdepthavg"_s]=sbal->z_s;
    (*face)["snowdepthavg_vert"_s]=sd_ver;

    sbal->input_rec1.S_n =sbal->input_rec2.S_n;
    sbal->input_rec1.I_lw =sbal->input_rec2.I_lw;
    sbal->input_rec1.T_a =sbal->input_rec2.T_a;
    sbal->input_rec1.e_a =sbal->input_rec2.e_a;
    sbal->input_rec1.u =sbal->input_rec2.u;
    sbal->input_rec1.T_g =sbal->input_rec2.T_g;
    sbal->input_rec1.ro =sbal->input_rec2.ro;

    // reset flag
//    g->dead = 0;
}

void snobal::checkpoint(mesh& domain,  netcdf& chkpt)
{

    chkpt.create_variable1D("snobal:m_s",domain->size_faces());
    chkpt.create_variable1D("snobal:rho",domain->size_faces());
    chkpt.create_variable1D("snobal:T_s",domain->size_faces());
    chkpt.create_variable1D("snobal:T_s_0",domain->size_faces());
    chkpt.create_variable1D("snobal:T_s_l",domain->size_faces());
    chkpt.create_variable1D("snobal:z_s",domain->size_faces());
    chkpt.create_variable1D("snobal:h2o_sat",domain->size_faces());
    chkpt.create_variable1D("snobal:max_h2o_vol",domain->size_faces());
    chkpt.create_variable1D("snobal:sum_runoff",domain->size_faces());
    chkpt.create_variable1D("snobal:sum_melt",domain->size_faces());
    chkpt.create_variable1D("snobal:E_s_sum",domain->size_faces());
    chkpt.create_variable1D("snobal:melt_sum",domain->size_faces());
    chkpt.create_variable1D("snobal:ro_pred_sum",domain->size_faces());
    chkpt.create_variable1D("snobal:h2o_total",domain->size_faces());

//netcdf puts are not threadsafe.
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        snodata *g = face->get_module_data<snodata>(ID);
        auto *sbal = &(g->data);

        chkpt.put_var1D("snobal:m_s",i,sbal->m_s);
        chkpt.put_var1D("snobal:rho",i,sbal->rho);
        chkpt.put_var1D("snobal:T_s",i,sbal->T_s);
        chkpt.put_var1D("snobal:T_s_0",i,sbal->T_s_0);
        chkpt.put_var1D("snobal:T_s_l",i,sbal->T_s_l);
        chkpt.put_var1D("snobal:z_s",i, sbal->z_s);
        chkpt.put_var1D("snobal:h2o_sat",i,sbal->h2o_sat);
        chkpt.put_var1D("snobal:max_h2o_vol",i,sbal->max_h2o_vol);

        chkpt.put_var1D("snobal:sum_runoff",i,g->sum_runoff);
        chkpt.put_var1D("snobal:sum_melt",i,g->sum_melt);
        chkpt.put_var1D("snobal:E_s_sum",i,sbal->E_s_sum);
        chkpt.put_var1D("snobal:melt_sum",i,sbal->melt_sum);
        chkpt.put_var1D("snobal:ro_pred_sum",i,sbal->ro_pred_sum);
        chkpt.put_var1D("snobal:h2o_total",i,sbal->h2o_total);
    }

}

void snobal::load_checkpoint(mesh& domain, netcdf& chkpt)
{
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        snodata *g = face->get_module_data<snodata>(ID);
        auto *sbal = &(g->data);

        sbal->m_s = chkpt.get_var1D("snobal:m_s",i);
        sbal->rho = chkpt.get_var1D("snobal:rho",i);
        sbal->T_s = chkpt.get_var1D("snobal:T_s",i);
        sbal->T_s_0 = chkpt.get_var1D("snobal:T_s_0",i);
        sbal->T_s_l = chkpt.get_var1D("snobal:T_s_l",i);
        sbal->z_s =  chkpt.get_var1D("snobal:z_s",i);
        sbal->h2o_sat =  chkpt.get_var1D("snobal:h2o_sat",i);
        sbal->max_h2o_vol = chkpt.get_var1D("snobal:max_h2o_vol",i);

        g->sum_runoff = chkpt.get_var1D("snobal:sum_runoff",i);
        g->sum_melt = chkpt.get_var1D("snobal:sum_melt",i);
        sbal->E_s_sum = chkpt.get_var1D("snobal:E_s_sum",i);
        sbal->melt_sum = chkpt.get_var1D("snobal:melt_sum",i);
        sbal->ro_pred_sum = chkpt.get_var1D("snobal:ro_pred_sum",i);
        sbal->h2o_total = chkpt.get_var1D("snobal:h2o_total",i);

        sbal->init_snow();
    }
}
