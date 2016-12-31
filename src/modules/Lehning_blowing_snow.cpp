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
    double aT = (1710.0+1.36*T)*pow(10.,-9.0);
    double Qt = aT * pow(u10,4);
    return Qt;
}

double F(double T)
{
    double D = 2.82*pow(10,-5); // diffusivity of water in air, m^2/s
//    mio::Atmosphere::
    double lamdba = 0.0243; // thermal conductivity of air, (W/(m K))
  //  PhysConst::Ls / ( lamdba + (T+273.)) * ((PhysConst::Ls * PhysConst::M) /(PhysConst::R*(T+273.)) - 1 ) + 1/(D*)
}
double Qs (double T, double u10, double RH)
{
    double sigma_2 = (1-RH)/100.;
    double b = 137.6/(25*25*25*25*25);

    return (b*sigma_2)/F(T)*pow(u10,5);
}

dbsm::dbsm(config_file cfg)
        :module_base(parallel::domain)
{
    depends("U_2m_above_srf");
    depends("vw_dir");
    depends("t");

    provides("dbsm_pq");
    provides("q_t");
    provides("sum_q_t");

}

void dbsm::init(mesh domain)
{
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        face->set_face_data("sum_q_t",0);
    }
}

void dbsm::run(mesh domain)
{
    arma::sp_mat A(domain->size_faces(), domain->size_faces());
    arma::vec b_subl(domain->size_faces(), arma::fill::zeros);
    arma::vec b(domain->size_faces(), arma::fill::zeros);

    double F = 1000.0;

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        double T = face->face_data("t");
        double phi = face->face_data("vw_dir");
        double u10 = face->face_data("U_2m_above_srf");
        //scale back up to u10
//        Atmosphere::log_scale_wind(U_R, Z_R, Z_2m_above_srf, snowdepthavg);

        double alpha_i = F/(6*face->get_area());

        arma::vec u(2);
        Vector_2 v = math::gis::bearing_to_cartesian(phi);
        u(0) = v.x(); //U_x
        u(1) = v.y(); //U_y

        arma::vec m1(2);
        m1(0) = face->edge_unit_normal(0).x();
        m1(1) = face->edge_unit_normal(0).y();


        arma::vec m2(2);
        m2(0) = face->edge_unit_normal(1).x();
        m2(1) = face->edge_unit_normal(1).y();


        arma::vec m3(2);
        m3(0) = face->edge_unit_normal(2).x();
        m3(1) = face->edge_unit_normal(2).y();


        //edge lengths
        double E1 = face->edge_length(0);
        double E2 = face->edge_length(1);
        double E3 = face->edge_length(2);

        A(i,face->cell_id) =  alpha_i*(E1*arma::dot(u,m1) + E2*arma::dot(u,m2)+E3*arma::dot(u,m3) +1/alpha_i );

        //0 flux BC
        auto neigh = face->neighbor(0);
        if(neigh != nullptr)
        {
            A(i, neigh->cell_id) = alpha_i * E1 * arma::dot(u, m1);
        }

        neigh = face->neighbor(1);
        if(neigh != nullptr)
        {
            A(i, neigh->cell_id) = alpha_i * E2 * arma::dot(u, m2);
        }

        neigh = face->neighbor(2);
        if(neigh != nullptr)
        {
            A(i, neigh->cell_id) = alpha_i * E3 * arma::dot(u, m3);
        }

        double pq = p10_dry(T, u10) * Qt(T, u10);
        b[i] = pq < 1e-6? 0 : -pq;

        face->set_face_data("dbsm_pq", b[i] );
    }


    arma::vec x = arma::spsolve(A,b);

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

        face->set_face_data("q_t", x(i)*global_param->dt());
        face->set_face_data("sum_q_t",  face->face_data("sum_q_t") + x(i)*global_param->dt());
    }
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
