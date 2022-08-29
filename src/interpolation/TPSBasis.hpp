#pragma once
#include <gsl/gsl_sf_expint.h>
#include <cmath>

// Needed if x>800 for old version. gsl_sf_expint_E1 overflows for x>800
//#include <gsl/gsl_errno.h>
//static gsl_error_handler_t *old_error_handler=gsl_set_error_handler_off();

template <typename T>
T TPSBasis(T x){
  T c = 0.577215;
  // gsl_sf_expint_E1(x) < epsilon_double and ln(x) > 1 for x > 32
  // so -(log(x) + c + gsl_sf_expint_E1(x)) == -(log(x) + c) is true for x > 32
  return (x <= 32) ? -(log(x) + c + gsl_sf_expint_E1(x)) : -(log(x) + c);
}

// TODO this function solely exists for timing purposes and should be removed
template <typename T>
T OLD_TPSBasis(T x){
  return -(log(x) + 0.577215 + gsl_sf_expint_E1(x));
}
