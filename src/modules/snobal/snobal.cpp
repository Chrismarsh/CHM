#include "snobal.hpp"

snobal::snobal(config_file cfg)
        : module_base(parallel::data)
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
    provides("dead");
    provides("iswr_net");
    provides("isothermal");
    provides("ilwr_out");
    provides("sum_snowpack_runoff");
    provides("sum_melt");
    provides("snowdepthavg");
}

void snobal::init(mesh domain)
{
    drift_density = cfg.get("drift_density",300.);

    //store all of snobals global_param variables from this timestep to be used as ICs for the next timestep
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
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
        sbal->relative_hts = 1;  //docs are wrong. 1 == absolute

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
            face->set_face_data("swe", sbal->m_s);
            face->set_face_data("R_n", sbal->R_n);
            face->set_face_data("H", sbal->H);
            face->set_face_data("E", sbal->L_v_E);
            face->set_face_data("G", sbal->G);
            face->set_face_data("M", sbal->M);
            face->set_face_data("dQ", sbal->delta_Q);
            face->set_face_data("cc", sbal->cc_s);
            face->set_face_data("T_s", sbal->T_s);
            face->set_face_data("T_s_0", sbal->T_s_0);
            face->set_face_data("iswr_net", sbal->S_n);
            face->set_face_data("isothermal", sbal->isothermal);
            face->set_face_data("ilwr_out", sbal->R_n - sbal->S_n - sbal->I_lw);
            face->set_face_data("snowmelt_int", 0);
            face->set_face_data("sum_melt", g->sum_melt);
            face->set_face_data("sum_snowpack_runoff", g->sum_runoff);

            face->set_face_data("snowdepthavg", sbal->z_s);
        }

        

    }
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

    double albedo = face->face_data("snow_albedo");

    // Optional inputs if there is a canopy or not
    double ilwr;
    if(has_optional("ilwr_subcanopy")) {
        ilwr = face->face_data("ilwr_subcanopy");
    } else {
        ilwr = face->face_data("ilwr");
    }

    // Optional inputs if there is a canopy or not
    double rh;
    if(has_optional("rh_subcanopy")) {
        rh = face->face_data("rh_subcanopy");
    } else {
        rh = face->face_data("rh");
    }

    // Optional inputs if there is a canopy or not
    double t;
    if(has_optional("ta_subcanopy")) {
        t = face->face_data("ta_subcanopy");
    } else {
        t = face->face_data("t");
    }

    double ea = mio::Atmosphere::vaporSaturationPressure(t+273.15)  * rh/100.;

    // Optional inputs if there is a canopy or not
    if(has_optional("iswr_subcanopy")) {
        sbal->input_rec2.S_n = (1-albedo) * face->face_data("iswr_subcanopy");
    } else {
        sbal->input_rec2.S_n = (1-albedo) * face->face_data("iswr");
    }
    sbal->input_rec2.I_lw = ilwr;
    sbal->input_rec2.T_a = t+FREEZE;
    sbal->input_rec2.e_a = ea;

    // Optional inputs if there is a canopy or not
    sbal->input_rec2.u = face->face_data("U_2m_above_srf");

    sbal->input_rec2.T_g = -4+FREEZE; //TODO: FIx this with a gflux estimate
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
        p = face->face_data("p_subcanopy");
    } else {
        p = face->face_data("p");
    }

    if(p > 0)
    {
        sbal->precip_now = 1;
        sbal->m_pp = p;

        // Optional inputs if there is a canopy or not
        if(has_optional("frac_precip_snow_subcanopy")) {
            sbal->percent_snow = face->face_data("frac_precip_snow_subcanopy");
        } else {
            sbal->percent_snow = face->face_data("frac_precip_snow");
        }
        sbal->rho_snow = 100.; //http://ccc.atmos.colostate.edu/pdfs/SnowDensity_BAMS.pdf
        sbal->T_pp = t+FREEZE;
        sbal->stop_no_snow=0;
    }
    else
    {
        sbal->precip_now = 0;
        sbal->stop_no_snow=0;
    }
  
    if(has_optional("drift_mass"))
    {

        double mass = face->face_data("drift_mass");
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
        g->delta_avalanche_snowdepth = face->face_data("delta_avalanche_snowdepth");
    }
    if(has_optional("delta_avalanche_mass")) {
        g->delta_avalanche_swe = face->face_data("delta_avalanche_mass");
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

    face->set_face_data("dead",g->dead);

    face->set_face_data("swe",sbal->m_s);

    face->set_face_data("R_n",sbal->R_n);
    face->set_face_data("H",sbal->H);
    face->set_face_data("E",sbal->L_v_E);
    face->set_face_data("G",sbal->G);
    face->set_face_data("M",sbal->M);
    face->set_face_data("dQ",sbal->delta_Q);
    face->set_face_data("cc",sbal->cc_s);
    face->set_face_data("T_s",sbal->T_s);
    face->set_face_data("T_s_0",sbal->T_s_0);
    face->set_face_data("iswr_net",sbal->S_n);
    face->set_face_data("isothermal",sbal->isothermal);
    face->set_face_data("ilwr_out", sbal->R_n - sbal->S_n - sbal->I_lw);
//    face->set_face_data("snowmelt_int",sbal->ro_predict);

    face->set_face_data("snowmelt_int",swe_diff);
    face->set_face_data("sum_melt",g->sum_melt);
    face->set_face_data("sum_snowpack_runoff",g->sum_runoff);

    face->set_face_data("snowdepthavg",sbal->z_s);

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

