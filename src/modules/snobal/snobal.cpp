#include "snobal.hpp"

snobal::snobal()
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


}

void snobal::init(mesh domain)
{
    //store all of snobals global variables from this timestep to be used as ICs for the next timestep
    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        snodata* g = face->make_module_data<snodata>(ID);
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
        sbal->max_h2o_vol = 0.0001;
        sbal->rho = 0.;
        sbal->T_s = -75. + FREEZE;
        sbal->T_s_0 = -75. + FREEZE; //assuming no snow
        sbal->T_s_l = -75. + FREEZE;
        sbal->z_s = 0.;

        sbal->KT_WETSAND = 0.08;

        sbal->ro_data = 0;

        sbal->max_z_s_0 = .1;
        sbal->h2o_total = 0;
        sbal->isothermal = 0;

        /// Heights
        sbal->z_0 = 0.001; //fix TODO: hardcode z_0
        sbal->z_T = 2.6;
        sbal->z_u = 2.96;
        sbal->z_g = 0.1;
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
        sbal->tstep_info[DATA_TSTEP].time_step = 3600; // 1hr TODO: fix hard code snobal timesteps
        sbal->tstep_info[DATA_TSTEP].intervals = 0;
        sbal->tstep_info[DATA_TSTEP].threshold = 60;
        sbal->tstep_info[DATA_TSTEP].output = 0;


        sbal->tstep_info[NORMAL_TSTEP].level = NORMAL_TSTEP;
        sbal->tstep_info[NORMAL_TSTEP].time_step = 3600;
        sbal->tstep_info[NORMAL_TSTEP].intervals = sbal->tstep_info[DATA_TSTEP].time_step /
                                                             sbal->tstep_info[NORMAL_TSTEP].time_step;;
        sbal->tstep_info[NORMAL_TSTEP].threshold = 60;
        sbal->tstep_info[NORMAL_TSTEP].output = 0;


        sbal->tstep_info[MEDIUM_TSTEP].level = MEDIUM_TSTEP;
        sbal->tstep_info[MEDIUM_TSTEP].time_step = 3600. / 4.; //15min
        sbal->tstep_info[MEDIUM_TSTEP].intervals = sbal->tstep_info[NORMAL_TSTEP].time_step /
                                                             sbal->tstep_info[MEDIUM_TSTEP].time_step;
        sbal->tstep_info[MEDIUM_TSTEP].threshold = 10;
        sbal->tstep_info[MEDIUM_TSTEP].output = 0;


        sbal->tstep_info[SMALL_TSTEP].level = SMALL_TSTEP;
        sbal->tstep_info[SMALL_TSTEP].time_step = 3600. / 60.;
        sbal->tstep_info[SMALL_TSTEP].intervals = sbal->tstep_info[MEDIUM_TSTEP].time_step /
                                                            sbal->tstep_info[SMALL_TSTEP].time_step;
        sbal->tstep_info[SMALL_TSTEP].threshold = 1;
        sbal->tstep_info[SMALL_TSTEP].output = 0;


        ////////

        sbal->init_snow();

        //////
        

    }
}
snobal::~snobal()
{

}

void snobal::run(mesh_elem &elem, boost::shared_ptr <global> global_param)
{
    //debugging
    auto id = elem->_debug_ID;

    bool run=false;


    auto hour = global_param->hour();
    auto day = global_param->day();
    auto month = global_param->month();


    //get the previous timesteps data out of the global store.
    snodata* g = elem->get_module_data<snodata>(ID);
    if(g->dead)
    {
        return;
    }
    auto* sbal = &(g->data);

    sbal->_debug_id = id;

    sbal->P_a = mio::Atmosphere::stdAirPressure( elem->get_z());

    double albedo = elem->face_data("snow_albedo");
    double ilwr = elem->face_data("ilwr");
    double rh = elem->face_data("rh");
    double t = elem->face_data("t");
    //double ea = mio::Atmosphere::saturatedVapourPressure(t+273.15) * rh/100.;
    double ea = mio::Atmosphere::waterSaturationPressure(t+273.15);

    sbal->input_rec2.S_n = (1-albedo) * elem->face_data("iswr"); //
    sbal->input_rec2.I_lw = ilwr;
    sbal->input_rec2.T_a = t+FREEZE;
    sbal->input_rec2.e_a = ea;
    sbal->input_rec2.u = elem->face_data("vw");
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


    double p = elem->face_data("p");

    if(p > 0)
    {
        sbal->precip_now = 1;
        sbal->m_pp = p;
        sbal->percent_snow = elem->face_data("frac_precip_snow");
        sbal->rho_snow = 100.; //http://ccc.atmos.colostate.edu/pdfs/SnowDensity_BAMS.pdf
        sbal->T_pp = t+FREEZE;
        sbal->stop_no_snow=0;
    }
    else
    {
        sbal->precip_now = 0;
        sbal->stop_no_snow=0;
    }

   // sbal->init_snow();
    try
    {
        if(! sbal->do_data_tstep() )
        {
            g->dead=1;
//            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("snobal died"));
        }
    }catch(...)
    {
        g->dead=1;
//        BOOST_THROW_EXCEPTION(module_error() << errstr_info ("snobal died"));
    }


    elem->set_face_data("dead",g->dead);

    elem->set_face_data("swe",sbal->m_s);

    elem->set_face_data("R_n",sbal->R_n);
    elem->set_face_data("H",sbal->H);
    elem->set_face_data("E",sbal->L_v_E);
    elem->set_face_data("G",sbal->G);
    elem->set_face_data("M",sbal->M);
    elem->set_face_data("dQ",sbal->delta_Q);
    elem->set_face_data("cc",sbal->cc_s);
    elem->set_face_data("T_s",sbal->T_s);
    elem->set_face_data("T_s_0",sbal->T_s_0);
    elem->set_face_data("iswr_net",sbal->S_n);
    elem->set_face_data("isothermal",sbal->isothermal);



    sbal->input_rec1.S_n =sbal->input_rec2.S_n;
    sbal->input_rec1.I_lw =sbal->input_rec2.I_lw;
    sbal->input_rec1.T_a =sbal->input_rec2.T_a;
    sbal->input_rec1.e_a =sbal->input_rec2.e_a;
    sbal->input_rec1.u =sbal->input_rec2.u;
    sbal->input_rec1.T_g =sbal->input_rec2.T_g;
    sbal->input_rec1.ro =sbal->input_rec2.ro;

}

void snobal::run(mesh domain, boost::shared_ptr <global> global_param)
{

}