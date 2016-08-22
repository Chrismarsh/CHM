#include "snobal.hpp"

snobal::snobal(config_file cfg)
        : module_base(parallel::data)
{
    depends("p_rain");
    depends("p_snow");
    depends("iswr");
    depends("rh");
    depends("t");
    depends("vw");
    depends("p");
    depends("ilwr");

    depends("snow_albedo");

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
    provides("snow_depth");


}

void snobal::init(mesh domain, boost::shared_ptr<global> global)
{
    //store all of snobals global variables from this timestep to be used as ICs for the next timestep
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        snodata* g = face->make_module_data<snodata>(ID);
        g->sum_runoff = 0;
        g->sum_melt = 0;
        auto* sbal = &(g->data);
        g->dead=0;
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
        sbal->tstep_info[DATA_TSTEP].time_step = global->dt(); //3600; // 1hr TODO: fix hard code snobal timesteps
        sbal->tstep_info[DATA_TSTEP].intervals = 0;
        sbal->tstep_info[DATA_TSTEP].threshold = 60;
        sbal->tstep_info[DATA_TSTEP].output = 0;


        sbal->tstep_info[NORMAL_TSTEP].level = NORMAL_TSTEP;
        sbal->tstep_info[NORMAL_TSTEP].time_step = global->dt();//3600;
        sbal->tstep_info[NORMAL_TSTEP].intervals = sbal->tstep_info[DATA_TSTEP].time_step /
                                                             sbal->tstep_info[NORMAL_TSTEP].time_step;;
        sbal->tstep_info[NORMAL_TSTEP].threshold = 60;
        sbal->tstep_info[NORMAL_TSTEP].output = 0;


        sbal->tstep_info[MEDIUM_TSTEP].level = MEDIUM_TSTEP;
        sbal->tstep_info[MEDIUM_TSTEP].time_step = global->dt()/4; //3600. / 4.; //15min
        sbal->tstep_info[MEDIUM_TSTEP].intervals = sbal->tstep_info[NORMAL_TSTEP].time_step /
                                                             sbal->tstep_info[MEDIUM_TSTEP].time_step;
        sbal->tstep_info[MEDIUM_TSTEP].threshold = 10;
        sbal->tstep_info[MEDIUM_TSTEP].output = 0;


        sbal->tstep_info[SMALL_TSTEP].level = SMALL_TSTEP;
        sbal->tstep_info[SMALL_TSTEP].time_step = global->dt()/ 60;//3600. / 60.;
        sbal->tstep_info[SMALL_TSTEP].intervals = sbal->tstep_info[MEDIUM_TSTEP].time_step /
                                                            sbal->tstep_info[SMALL_TSTEP].time_step;
        sbal->tstep_info[SMALL_TSTEP].threshold = 1;
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

        //////
        

    }
}
snobal::~snobal()
{

}

void snobal::run(mesh_elem &face, boost::shared_ptr<global> global_param)
{
    //debugging
    auto id = face->_debug_ID;

    bool run=false;


    auto hour = global_param->hour();
    auto day = global_param->day();
    auto month = global_param->month();


    //get the previous timesteps data out of the global store.
    snodata* g = face->get_module_data<snodata>(ID);
    auto* sbal = &(g->data);

    sbal->_debug_id = id;

    sbal->P_a = mio::Atmosphere::stdAirPressure( face->get_z());

    double albedo = face->face_data("snow_albedo");
    double ilwr = face->face_data("ilwr");
    double rh = face->face_data("rh");
    double t = face->face_data("t");
    double ea = mio::Atmosphere::waterSaturationPressure(t+273.15)  * rh/100.;

    sbal->input_rec2.S_n = (1-albedo) * face->face_data("iswr");
    sbal->input_rec2.I_lw = ilwr;
    sbal->input_rec2.T_a = t+FREEZE;
    sbal->input_rec2.e_a = ea;
    sbal->input_rec2.u = face->face_data("vw");
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


    double p = face->face_data("p");

    if(p > 0)
    {
        sbal->precip_now = 1;
        sbal->m_pp = p;
        sbal->percent_snow = face->face_data("frac_precip_snow");
        sbal->rho_snow = 100.; //http://ccc.atmos.colostate.edu/pdfs/SnowDensity_BAMS.pdf
        sbal->T_pp = t+FREEZE;
        sbal->stop_no_snow=0;
    }
    else
    {
        sbal->precip_now = 0;
        sbal->stop_no_snow=0;
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
        auto details = "("+std::to_string(face->center().x()) + "," + std::to_string(face->center().y())+","+std::to_string(face->center().z())+")";
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

    face->set_face_data("snow_depth",sbal->z_s);

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

void snobal::run(mesh domain, boost::shared_ptr <global> global_param)
{

}