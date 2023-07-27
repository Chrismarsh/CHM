#pragma once
#include <boost/math/special_functions/expint.hpp>
#include <gsl/gsl_sf_expint.h>
#include <cmath>

/* Basis function used for building a thin plate spline. boost::math::expint is used here instead
 * of gsl_sf_expint because it automatically works with Boost's automatic differentiation */
template <typename T>
T TPSBasis(T x){
  T gamma = 0.5772156649015328606;
  // gsl_sf_expint_E1(x) < epsilon_double and ln(x) > 1 for x > 32
  // so -(log(x) + c + gsl_sf_expint_E1(x)) == -(log(x) + c) is true for x > 32
  return (x <= 32) ? -(log(x) + gamma + boost::math::expint(1,x)) : -(log(x) + gamma);
}


/* TODO OLD_TPSBasis exists solely for timing FunC's LUTs against direct evaluation. It can be safely removed.
 * The following two lines are needed if x>800 for OLD_TPSBasis because gsl_sf_expint_E1 overflows for x>800;
 * however, old_error_handler is set to gsl_set_error_handler_off() somewhere else in CHM */
//#include <gsl/gsl_errno.h>
//static gsl_error_handler_t *old_error_handler=gsl_set_error_handler_off();
template <typename T>
T OLD_TPSBasis(T x){
  return -(log(x) + 0.577215 + gsl_sf_expint_E1(x));
}
