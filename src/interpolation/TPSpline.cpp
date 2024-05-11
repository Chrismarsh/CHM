//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//


#include "TPSpline.hpp"

//#include <iostream>
//#define FUNC_DEBUG
#include <func/func.hpp>
#include "TPSBasis.hpp"
#include <boost/predef.h>

// Build FunC lookup table for -(log(x)+gamma+gsl_sf_expint_E1(x))

// hardcoded error of 1e-8
// There is a compiler bug fixed in gcc that interfers with boost 1.85.0 math and frounding-math (needed for cgal)
// so use the uniform spaced LUT to avoid the proplematic call
// bug is fixed in  12.4/13.3/14.0
// see https://github.com/boostorg/math/issues/1133 and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=109359
#if BOOST_COMP_GNUC && BOOST_COMP_GNUC < BOOST_VERSION_NUMBER(13,3,0)
static func::FailureProofTable<func::UniformEqSpaceInterpTable<3,double>,double> TPSBasis_LUT({TPSBasis<double>}, {1e-5, 32, 0.11111});
#else
static func::FailureProofTable<func::NonUniformEqSpaceInterpTable<3,double>,double> TPSBasis_LUT({FUNC_SET_F(TPSBasis,double)}, {1e-5, 32, 0.154207});
#endif

//static func::DirectEvaluation<double> TPSBasis_LUT({OLD_TPSBasis<double>}, 1e-8, 6500);


/* Each MPI process will print that text from the destructor after it finishes */
//static struct Counter{
//  double total = 0.0;
//  ~Counter(){ std::cerr << "Spent " << total << " seconds in spline::operator()" << std::endl; }
//} counter;


double thin_plate_spline::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)
{
    //func::Timer t;
    //see if we can reuse our
    if(sample_points.size() + 1 != size)
    {
        size = sample_points.size();
        size++; // need to make room for the physics
        A = MatrixXXd::Zero(size,size);
        b = VectorXd::Zero(size);
        x = VectorXd::Zero(size);
    }

    if(uninit_lu_decomp)
    {
        //build the LU decomp
        for (unsigned int i = 0; i < size - 1; i++)
        {
            double sxi = sample_points.at(i).get<0>(); //x
            double syi = sample_points.at(i).get<1>(); //y

            for (unsigned int j = i; j < size - 1; j++)
            {
                double sxj = sample_points.at(j).get<0>(); //x
                double syj = sample_points.at(j).get<1>(); //y

                double xdiff = (sxi - sxj);
                double ydiff = (syi - syj);

                //don't add in a duplicate point, otherwise we get nan
                if (xdiff == 0. && ydiff == 0.)
                    continue;

                double Rd = 0.;
                if (j == i) // diagonal
                {
                    Rd = 0.0;
                } else
                {
                    double dij = sqrt(xdiff * xdiff + ydiff * ydiff); //distance between this set of observation points

                    //none of the books and papers, despite citing Helena Mitášová, Lubos Mitáš seem to agree on the exact formula
                    //so I am following http://link.springer.com/article/10.1007/BF00893171#page-1
                    // eqn 10

                    dij = (dij * weight / 2.0) * (dij * weight / 2.0);

                    //Chang 4th edition 2008 uses bessel_k0
                    //gsl_sf_bessel_K0
                    // and has a -0.5 weight out fron
//                     Rd = -0.5/(pi*weight*weight)*( log(dij*weight/2.0) + c + gsl_sf_bessel_K0(dij*weight));

                    //And Hengl and Evans in geomorphometry p.52 do not, but have some undefined omega_0/omega_1 weights
                    //it is all rather confusing. But this follows Mitášová exactly, and produces essentially the same answer
                    //as the worked example in box 16.2 in Chang
                    // set Rd = -(log(dij) + c + gsl_sf_expint_E1(dij))
                    Rd = TPSBasis_LUT(dij);
                }

                A(i, j + 1) = Rd;
                A(j, i + 1) = Rd;

            }
        }


        //set physics and build b values
        for (unsigned int i = 0; i < size; i++)
        {
            A(i, 0) = 1;
            A(size - 1, i) = 1;
        }
        A(size - 1, 0) = 0;
        lu.compute(A);
    }


    //this will skip the above computation of the lu decomp next time this is called
    if(reuse_LU)
        uninit_lu_decomp = false;

    for(size_t i=0;i<size-1;i++)
    {
        b(i) = sample_points.at(i).get<2>() ;
    }

    b(size-1) = 0.0; //constant

    //solve equation
    int s				= 0;

    x =  lu.solve(b) ; //ldlt.solve(b);


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
        // set Rd equal to -(log(dij) + c + gsl_sf_expint_E1(dij))
         
        double Rd = TPSBasis_LUT(dij);

        z0 = z0 + x(i)*Rd;
    }

    //t.stop();
    //counter.total += t.duration();
    return z0;
}

thin_plate_spline::thin_plate_spline(size_t sz, std::map<std::string,std::string> config )
: thin_plate_spline()
{
    size = sz;
    size++; // need to make room for the physics

    A = MatrixXXd::Zero(size,size);
    b = VectorXd::Zero(size);
    x = VectorXd::Zero(size);

    auto itr = config.find("reuse_LU");
    if(itr != config.end())
    {
        if(config["reuse_LU"] == "true")
            reuse_LU = true;
    }
    reuse_LU = false;
}

thin_plate_spline::thin_plate_spline()
{
    pi          = 3.14159;
    c           = 0.577215; //euler constant
    weight      = 0.01;
    size        = 0;

    reuse_LU    = false;
    uninit_lu_decomp = true;
}

thin_plate_spline::~thin_plate_spline(){}
