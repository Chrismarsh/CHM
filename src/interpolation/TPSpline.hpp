#pragma once
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_expint.h>

#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>
#include <armadillo>

#include <boost/throw_exception.hpp>
#include <exception.hpp>
#include <iostream>
#include "interp_base.hpp"
#include "logger.hpp"

/**
* \class thin_plate_spline
*
* Thin plate spline with tensions interpolation
*/
class thin_plate_spline : public interp_base
{
public:
    thin_plate_spline();
    ~thin_plate_spline();

    /**
     * Allocating A,b,x, and p every time this is called adds up.
     * So the sz can be pre-set to pre allocate all the arrays.
     */
    thin_plate_spline(size_t sz,std::map<std::string,std::string> config = std::map<std::string,std::string>());

    /**
    * Spline of the sample_points at the query_point location.
    * \param sample_points Tuple of x,y,z values that comprise the sample points from which to interpolate from
    * \param query_point Tuple of x,y,z value that is the point to interpolate to
    * \return Interpolated value at the query_point
    */
    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point);

    bool reuse_LU;
private:
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> VectorXd;
    typedef Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> MatrixXXd;

    MatrixXXd A ;
    VectorXd b; // known values - constant value of 0 goes in b[size-1]
    VectorXd x;

    Eigen::FullPivLU< Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> > lu;
    double pi;
    double c; //euler constant
    double weight;
    bool uninit_lu_decomp;
    size_t size;

};
