#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */

/* Problem Constants */

#define NEQ   2                /* number of equations  */
#define Y1    RCONST(0.75)      /* initial y components */
#define Y2    RCONST(0.5)


#define RTOL  RCONST(1.0e-6)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-4)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-4)


#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(60.0)      /* first output time      */

#define TMULT RCONST(1.0)     /* output time factor     */
#define NOUT  1          /* number of output times */

struct module_bucket_face :  face_info
{
    module_bucket_face()
    {
        IC = new double[2];
        IC[0]=0;
        IC[1]=0;
    }

    double* IC;
};

/**
* \addtogroup modules
* @{
* \class leaky bucket
* \brief Calculates 2-bucket leaky bucket
*
* Calculates 2-buck leaky bucket
*
* Depends:
* -
*
* Provides:
* -
*/
class leaky_bucket : public module_base
{
public:
    leaky_bucket(std::string ID);
    ~leaky_bucket();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

/**
@}
*/