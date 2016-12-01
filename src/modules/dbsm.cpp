#include "dbsm.hpp" 


double p10_dry(double T, double u10)
{
    //eqn 10
    double snowage = 1; //hours
    double u_bar = 11.2 + 0.356*T+0.00706*T*T+0.8*log(1+snowage);

    //eqn 11
    double var = 4.3 + 0.145*T+0.00196*T*T;

    //eqn 12
    double p = 1/(1+exp( sqrt(M_PI) * (u_bar - u10) / var));

    return p;
}

double Qt (double T, double u10)
{
    double aT = (1710+1.36*T)*pow(10.,-9);
    double Qt = aT * pow(u10,4);
    return Qt;
}

//int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
//{
//    //magnitude of the flux
//    realtype qi;
//    //IC of magnitude of flux
//    qi = Ith(y,1);
//
//    dbsm::cvode_data* d = (dbsm::cvode_data*) user_data;
//
//    double F = 1000.;
//    double u_bar = 11.2 + 0.365*T + 0.00196*T*T;
//    // x
//    double sum_x =
//
//    Ith(ydot,1) = p10_dry(T,U10) * Qt(T,U10) - F/(3*Vi)*sum_x;
//
//    return(0);
//}



dbsm::dbsm(config_file cfg)
        :module_base(parallel::data)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("t");

    provides("dbsm_u10");
    provides("q_t");

}

void dbsm::init(mesh domain)
{

//#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<face_data>(ID);

    }
}

void dbsm::run(mesh_elem& face)
{
    double T = face->face_data("t");
    double phi = face->face_data("vw_dir");
    double u10 = 10*face->face_data("U_2m_above_srf");

    arma::vec n(2);
    n(0)  = cos(phi * 3.14159/180.0); //U_x
    n(1)  = -sin(phi * 3.14159/180.0); //U_y

    double F = 1000.0;

    auto d  = face->get_module_data<face_data>(ID);

    //current q_x, q_y components to use as ICs
    arma::vec q(2);
    q = d->IC;

    //magnitude of previous time step's value, transport
    double q_t_prev = arma::norm(q);

    arma::vec result(2);
    result(0)=0;
    result(1)=0;

    for(int i = 0; i <3; i++)
    {
        arma::vec q_xn(2);
        auto neigh = face->neighbor(i);
        double E = 100; //TODO: fix
        if(neigh != nullptr)
        {
            //make unit wind vector at our neighbour
            auto dir = neigh->face_data("vw_dir");
            //make vector of trans flux
            q_xn = n * q_t_prev; //elem wise mult


            result += E*(q_xn - q) / math::gis::distance(neigh->center(),face->center());
        }
        else
        {
            q_xn(0) = 0;
            q_xn(1) = 0;

//            result += E*(q_xn - q) / math::gis::distance(neigh->center(),face->center());
        }

    }

    double scaled = p10_dry(T,u10)*Qt(T,u10);
    double trans = F/(3.*face->get_area()) * arma::dot(n,result);

    double blowing = (p10_dry(T,u10)*Qt(T,u10) - F/(3*face->get_area()) * arma::dot(n,result)) * global_param->dt() + q_t_prev;

   d->IC = n*blowing;

    face->set_face_data("dbsm_u10",u10);
    face->set_face_data("q_t",blowing);
}

//void dbsm::run(mesh_elem& face)
//{
//
//    realtype reltol, t, tout;
//    N_Vector y, abstol;
//    void *cvode_mem;
//    int flag, flagr, iout;
//
//    cvode_data* d = new cvode_data;
//
//    double phi = face->face_data("vw_dir");
//
//    d->T = face->face_data("t");
//    d->u10 = face->face_data("2m_sufrace_wind");
//
//    d->n(0)  = cos(phi * 3.14159/180.0); //U_x
//    d->n(1)  = -sin(phi * 3.14159/180.0); //U_y
//
//
//    userdata[3] = face->get_area(); //area
//    userdata[4] = 1;//edge length
//    userdata[5] = 1;//edge length
//    userdata[6] = 1;//edge length
//
//    auto n1 = face->neighbor(0);
//    auto n2 = face->neighbor(1);
//    auto n3 = face->neighbor(2);
//
//    // assume Dirichlet = 0 BC
//    userdata[7] = 0;
//    userdata[8] = 0;
//
//    if(n1)
//        userdata[7] = face->face_data("q_x");
//    if(n2)
//        userdata[8] = face->face_data("q_x");
//
//    y = abstol = NULL;
//    cvode_mem = NULL;
//
//    /* Create serial vector of length NEQ for I.C. and abstol */
//    y = N_VNew_Serial(NEQ);
//
//    abstol = N_VNew_Serial(NEQ);
//
//
//    module_bucket_face* info = nullptr;
//
//
//    face_info* tmp_info = face->module_face_data(ID);
//
//
//        info->IC[0]=Y1;
//        info->IC[1]=Y2;
//
//
//    /* Initialize y */
//    Ith(y,1) = info->IC[0];
//    Ith(y,2) = info->IC[1];
//
//    /* Set the scalar relative tolerance */
//    reltol = RTOL;
//    /* Set the vector absolute tolerance */
//    Ith(abstol,1) = ATOL1;
//    Ith(abstol,2) = ATOL2;
//
//    /* Call CVodeCreate to create the solver memory and specify the
//    * Backward Differentiation Formula and the use of a Newton iteration */
//    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
//    /* Call CVodeInit to initialize the integrator memory and specify the
//     * user's right hand side function in y'=f(t,y), the inital time T0, and
//     * the initial dependent variable vector y. */
//    flag = CVodeInit(cvode_mem, fn, T0, y);
//
//
//    /* Call CVodeSVtolerances to specify the scalar relative tolerance
//     * and vector absolute tolerances */
//    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
//
//
//    /* Call CVDense to specify the CVDENSE dense linear solver */
//    flag = CVDense(cvode_mem, NEQ);
//
//    flag = CVodeSetUserData(cvode_mem,userdata);
//
//    flag = CVodeSetMaxNumSteps(cvode_mem,9001);
//
//    /* In loop, call CVode, print results, and test for error.
//    Break out of loop when NOUT preset output times have been reached.  */
//    iout = 0;
//    tout = T1;
////    while(1) {
//    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
////        LOG_DEBUG << "t="<<t<<" Ith="<<Ith(y,1); //<<" Ith2=" << Ith(y,2);
//
////        check_flag(&flag, "CVode", 1)
////        if (flag == CV_SUCCESS) {
////            iout++;
////            tout *= TMULT;
////        }
////
////        if (iout == NOUT) break;
////    }
//
//    elem->set_face_data("z1",Ith(y,1));
//    elem->set_face_data("z2",Ith(y,2));
//
//    info->IC[0] = Ith(y,1);
//    info->IC[1] = Ith(y,2);
//    elem->set_module_face_data(ID, (info));
//
//    /* Free y and abstol vectors */
//    N_VDestroy_Serial(y);
//    N_VDestroy_Serial(abstol);
//
//    /* Free integrator memory */
//    CVodeFree(&cvode_mem);
//
//    delete[] userdata;
//}

dbsm::~dbsm()
{

}
