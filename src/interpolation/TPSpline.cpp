
#include "TPSpline.hpp"

double thin_plate_spline::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)
{
    if(!static_size)
    {
        size = sample_points.size();
        size++; // need to make room for the constants
        A.resize(size,size);
        b.resize(size);

    }


    for(unsigned int i=0;i < size-1; i++)
    {
        double sxi = sample_points.at(i).get<0>(); //x
        double syi = sample_points.at(i).get<1>(); //y

        for(unsigned int j=i; j < size-1; j++)
        {
            double sxj = sample_points.at(j).get<0>(); //x
            double syj = sample_points.at(j).get<1>(); //y

            double xdiff = (sxi - sxj);
            double ydiff = (syi - syj);

            //don't add in a duplicate point, otherwise we get nan
            if(xdiff == 0. && ydiff == 0.)
                continue;

            double Rd=0.;
            if (j == i) // diagonal
            {
                Rd = 0.0;
            } else
            {
                double dij = sqrt(xdiff*xdiff + ydiff*ydiff); //distance between this set of observation points

                //none of the books and papers, despite citing Helena Mitášová, Lubos Mitáš seem to agree on the exact formula
                //so I am following http://link.springer.com/article/10.1007/BF00893171#page-1

                dij = (dij * weight/2.0) * (dij * weight/2.0);

                //Chang 4th edition 2008 uses bessel_k0
                //gsl_sf_bessel_K0
                // and has a -0.5 weight out fron
                // Rd = -0.5/(pi*weight*weight)*( log(dij*weight/2.0) + c + gsl_sf_bessel_K0(dij*weight));

                //And Hengl and Evans in geomorphometry p.52 do not, but have some undefined omega_0/omega_1 weights
                //it is all rather confusing. But this follows Mitášová exactly, and produces essentially the same answer
                //as the worked example in box 16.2 in Chang
                Rd = -(log(dij) + c + gsl_sf_expint_E1(dij));

            }

            A(i,j+1) = Rd;
            A(j,i+1) = Rd;

        }
    }

    //set constants and build b values
    for(unsigned int i = 0; i< size;i++)
    {
        A(i,0) = 1;
//        A[i*size+0] = 1;

        A(size-1,i)=1;
//        A[(size-1)*size+i] = 1;
    }
    A(size-1,0)=0;
//    A[(size-1)*(size)+0] = 0.0;

    for(size_t i=0;i<size-1;i++)
    {
        b(i) = sample_points.at(i).get<2>() ;
    }
//        B[i] = sample_points.at(i).get<2>() ;

    b(size-1) = 0.0; //constant

    //solve equation
    int s				= 0;

    arma::vec x = arma::solve(A,b);
    double z0 = x(0);//little a

    //skip x[0] we already pulled off above
    double ex = query_point.get<0>();
    double ey =  query_point.get<1>();

    for (unsigned int i = 1; i < x.size() ;i++)
    {
        double sx = sample_points.at(i-1).get<0>(); //x
        double sy = sample_points.at(i-1).get<1>(); //y

        double xdiff = (sx  - ex);
        double ydiff = (sy  - ey);
        double dij = sqrt(xdiff*xdiff + ydiff*ydiff);
        dij = (dij * weight/2.0) * (dij * weight/2.0);
        double Rd = -(log(dij) + c + gsl_sf_expint_E1(dij));

        z0 = z0 + x(i)*Rd;
    }

    return z0;
}

thin_plate_spline::thin_plate_spline(size_t sz)
: thin_plate_spline()
{

    size = sz;
    size++; // need to make room for the constants

    A.resize(size,size );
    b.resize(size); // known values - constant value of 0 goes in b[size-1]

    static_size = true;


}
thin_plate_spline::thin_plate_spline()
{
    pi          = 3.14159;
    c           = 0.577215; //euler constant
    weight      = 0.1;
    size        = 0;
    static_size = false;
}
thin_plate_spline::~thin_plate_spline()
{

}