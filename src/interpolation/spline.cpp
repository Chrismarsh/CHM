#include "spline.hpp"

double spline::operator()(station_list&  stations, mesh_elem& elem,boost::shared_ptr<interp_visitor> visitor, boost::shared_ptr<global> global_param)
{

    unsigned int size = stations.size();
    size++; // need to make room for the constants

    double *A = new double[ size * size ];
    double pi = 3.14159;
    double c = 0.577215; //euler constant
    double *B = new double[size]; // known values - constant value of 0 goes in b[size-1]
    double weight = 0.1;
    for(unsigned int i=0;i < size-1; i++)
    {
            for(unsigned int j=i; j < size-1; j++)
            {
                double sxi = stations.at(i)->get_x();
                double syi = stations.at(i)->get_y();
                
                double sxj = stations.at(j)->get_x();
                double syj = stations.at(j)->get_y();

                double xdiff = (sxi  - sxj);
                double ydiff = (syi  - syj);
                double dij = pow(sqrt(
                                     pow(xdiff,2.0) + pow(ydiff,2.0)
                                    )
                                ,2.0
                            );
                //distance between all the observation points
                double Rd = 0.0;
               
                if ( (dij - 0.0) < 0.000001)
                        Rd = 0.0;
                else
                        Rd = -1.0/(2.0*pi*pow(weight,2.0))*( log(dij*weight/2.0) + c + gsl_sf_bessel_K0(dij*weight));

                //A(i,j+1) = Rd;
                A[i * size + (j+1)] = Rd;

                //A(j,i+1) = Rd;
                A[j * size + (i+1)] = Rd;			
            }
    }

    //set constants and build b values
    for(unsigned int i = 0; i< size;i++)
    {
            A[i*size+0] = 1;
            A[(size-1)*size+i] = 1;
    }

    int i=0;
    for(auto& itr : stations)
    {
            B[i] = visitor->lower(elem,itr,global_param);
            ++i;
    }
    B[size-1] = 0.0; //constant


    A[(size-1)*(size)] = 0.0;

    //solve equation via gsl via lu decomp
    gsl_matrix_view m 	= gsl_matrix_view_array (A, size, size);
    gsl_vector_view b	= gsl_vector_view_array (B, size);
    gsl_vector *x		= gsl_vector_alloc (size);
    int s				= 0;
    gsl_permutation * p = gsl_permutation_alloc (size);
    gsl_linalg_LU_decomp (&m.matrix, p, &s);
    gsl_linalg_LU_solve (&m.matrix, p, &b.vector, x);

    //location of the unknown
    unsigned int ux = elem.get_x();
    unsigned int uy = elem.get_y();

    double z0 = x->data[0];//little a

    //-1 to skip the x[0] we already pulled off above
    for (unsigned int i = 0;i < (x->size - 1) ;i++)
    {
        double sx = stations.at(i)->get_x();
        double sy = stations.at(i)->get_y();
        double ex = elem.get_x();
        double ey = elem.get_y();
        double xdiff = (sx  - ex);
        double ydiff = (sy  - ey);
        double d = pow(sqrt(
                             pow(xdiff,2.0) + pow(ydiff,2.0)
                            )
                        ,2.0
                    );
            double Rd = 0.0;

            if ((d - 0.0) < 0.000001)
                    Rd = 0.0;
            else
                    Rd = -1.0/(2.0*pi*pow(weight,2.0))*( log(d*weight/2.0) + c + gsl_sf_bessel_K0(d*weight));

            z0 = z0 + x->data[i+1]*Rd; 
    }


    gsl_permutation_free (p);
    gsl_vector_free (x);
    delete B;
    delete A;

    z0 = visitor->raise(z0,elem,global_param);
    return z0;
}
spline::spline()
{
    
}

spline::~spline()
{
    
}