#pragma once
//#include <boost/math/special_functions/expint.hpp>
//#include <boost/math/differentiation/autodiff.hpp>
#include <gsl/gsl_sf_expint.h>
#include <cmath>
#define FUNC(x) -(log(x) + 0.577215 + gsl_sf_expint_E1(x))
#define FUNCNAME "CHM Elliptic Integral Equation"

//static std::uniform_real_distribution<double> my_noise(-1,1);
//static std::default_random_engine my_gen;

template <typename T>
T CHM_Elliptic_Equation(T x){
  return -(log(x) + 0.577215 + gsl_sf_expint_E1(x));
  //return -(log(x) + 0.577215 + boost::math::expint(1,x));
  //return -(log(x) + 0.577215 + boost::math::expint(1,x))*(1.0+1e-8*my_noise(my_gen));
}
