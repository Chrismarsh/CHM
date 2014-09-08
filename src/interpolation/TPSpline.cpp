
#include "TPSpline.hpp"

double thin_plate_spline::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_points )
{

    size_t size = sample_points.size();
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
            double sxi = sample_points.at(i).get<0>(); //x
            double syi = sample_points.at(i).get<1>(); //y

            double sxj = sample_points.at(j).get<0>(); //x
            double syj = sample_points.at(j).get<1>(); //y

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

//    int i=0;
//    for(auto& itr : stations)
//    {
//        B[i] = visitor->lower(elem,itr,global_param);
//        ++i;
//    }

    for(size_t i=0;i<size-1;i++)
        B[i] = sample_points.at(i).get<2>() ; //z

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

    double z0 = x->data[0];//little a

    //-1 to skip the x[0] we already pulled off above
    for (unsigned int i = 0;i < (x->size - 1) ;i++)
    {
        double sx = sample_points.at(i).get<0>();
        double sy = sample_points.at(i).get<1>();
        double ex = query_points.get<0>();
        double ey =  query_points.get<1>();
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
    delete[] B;
    delete[] A;

//    z0 = visitor->raise(z0,elem,global_param);
    return z0;


}

thin_plate_spline::thin_plate_spline()
{

}
thin_plate_spline::~thin_plate_spline()
{

}