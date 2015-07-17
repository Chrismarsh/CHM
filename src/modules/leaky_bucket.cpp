#include "leaky_bucket.hpp"

static int check_flag(void *flagvalue, char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

        /* Check if flag < 0 */
    else if (opt == 1) {
        errflag = (int *) flagvalue;
        if (*errflag < 0) {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                    funcname, *errflag);
            return(1); }}

        /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
        fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                funcname);
        return(1); }

    return(0);
}

static int fn(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
    realtype y1, y2, ydot1;

    y1 = Ith(y,1);
    y2 = Ith(y,2);

    double* data = (double*) user_data;

    double P = data[0];
    double PET = data[1];
    double k = data[2];
    double k2= data[3];
    double LAI = data[4];
    double f = data[5];

    double z1 = y1;
    double z2 = y2;

    //non-leaky
    Ith(ydot,1) = P-PET*0.01*((5.0*z1-2.0*pow(z1,2.0))*(1.0/3.0))-P* pow(z1,(1.0/2.0)*LAI);
    if ( (z1-1.0) < 0.001)
    {
        //redo it and couple the eqns
        Ith(ydot, 1) = P - PET * 0.01 * ((5.0 * z1 - 2.0 * pow(z1, 2.0)) * (1.0 / 3.0)) - P * pow(z1, (1.0 / 2.0) * LAI) - (1.0 - f) * k * pow(z1, 2.0);
        Ith(ydot, 2) = (1.0 - f) * k * pow(z1, 2.0) - k2 * pow(z2, 2.0);
    }


//    Ith(ydot,1) = P-PET*0.01*((5.0*z1-2.0*pow(z1,2.0))*(1.0/3.0))-P* pow(z1,(1.0/2.0)*LAI)-f*k-(1.0-f)*k*pow(z1,2.0);
//    Ith(ydot,2) = (1.0-f)*k*pow(z1,2.0)-k2*pow(z2,2.0);



    return(0);
}



leaky_bucket::leaky_bucket( std::string ID)
{

    _depends->push_back("p");
    _depends->push_back("ET");


    _provides->push_back("z1");
    _provides->push_back("z2");

    this->ID = ID;
    _parallel_type = parallel::data;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

void leaky_bucket::run(mesh_elem& elem, boost::shared_ptr<global> global_param)
{



    realtype reltol, t, tout;
    N_Vector y, abstol;
    void *cvode_mem;
    int flag, flagr, iout;

    double* userdata = new double[6];

    userdata[0] = elem->face_data("p");
    userdata[1] = elem->face_data("ET");
    userdata[2] = 0.000004650819; //k1 mm/hr 0.0021
    userdata[3] = 0.00000000007067933; //k2 mm/h unsat sand
    userdata[4] = 2.0; //lai
    userdata[5] = 0.4; //f

    y = abstol = NULL;
    cvode_mem = NULL;

    /* Create serial vector of length NEQ for I.C. and abstol */
    y = N_VNew_Serial(NEQ);

    abstol = N_VNew_Serial(NEQ);


    module_bucket_face* info = nullptr;


    face_info* tmp_info = elem->module_face_data(ID);

    if (tmp_info == nullptr) //need to create
    {
        info = new module_bucket_face();
        info->IC[0]=Y1;
        info->IC[1]=Y2;
    }
    else
    {
        info = reinterpret_cast<module_bucket_face*>(tmp_info);
    }


    /* Initialize y */
    Ith(y,1) = info->IC[0];
    Ith(y,2) = info->IC[1];

    /* Set the scalar relative tolerance */
    reltol = RTOL;
    /* Set the vector absolute tolerance */
    Ith(abstol,1) = ATOL1;
    Ith(abstol,2) = ATOL2;

    /* Call CVodeCreate to create the solver memory and specify the
    * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function in y'=f(t,y), the inital time T0, and
     * the initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, fn, T0, y);


    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol, abstol);


    /* Call CVDense to specify the CVDENSE dense linear solver */
    flag = CVDense(cvode_mem, NEQ);

    flag = CVodeSetUserData(cvode_mem,userdata);

    flag = CVodeSetMaxNumSteps(cvode_mem,9001);

    /* In loop, call CVode, print results, and test for error.
    Break out of loop when NOUT preset output times have been reached.  */
    iout = 0;
    tout = T1;
//    while(1) {
        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
//        LOG_DEBUG << "t="<<t<<" Ith="<<Ith(y,1); //<<" Ith2=" << Ith(y,2);

//        check_flag(&flag, "CVode", 1)
//        if (flag == CV_SUCCESS) {
//            iout++;
//            tout *= TMULT;
//        }
//
//        if (iout == NOUT) break;
//    }

    elem->set_face_data("z1",Ith(y,1));
    elem->set_face_data("z2",Ith(y,2));

    info->IC[0] = Ith(y,1);
    info->IC[1] = Ith(y,2);
    elem->set_module_face_data(ID, (info));

    /* Free y and abstol vectors */
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(abstol);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    delete[] userdata;
}

leaky_bucket::~leaky_bucket()
{

}