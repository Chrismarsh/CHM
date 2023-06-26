#pragma once
#include <gsl/gsl_sf_expint.h>
#include <cmath>

// the following two lines disables the overflow exception gsl_sf_expint_E1 throws for x>800
//#include <gsl/gsl_errno.h>
//static gsl_error_handler_t *old_error_handler=gsl_set_error_handler_off();

template <typename T>
T TPSBasis(T x){
  T gamma = 0.57721566490153286060651209008240243104215933593992;
  // gsl_sf_expint_E1(x) < epsilon_double and ln(x) > 1 for x > 32
  // so -(log(x) + gamma + gsl_sf_expint_E1(x)) == -(log(x) + gamma) is true for x > 32
  return (x <= 32) ? -(log(x) + gamma + gsl_sf_expint_E1(x)) : -(log(x) + gamma);
}

// this function only exists for profiling
template <typename T>
T OLD_TPSBasis(T x){
  return -(log(x) + 0.577215 + gsl_sf_expint_E1(x));
}
