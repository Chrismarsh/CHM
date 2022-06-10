#pragma once
#include <boost/math/differentiation/autodiff.hpp>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_expint.h> // need libgsl (gnu scientific library)
#include <cmath>
#include <iostream>
#define FUNC(x) -(log(x) + 0.577215 + gsl_sf_expint_E1(x))
#define FUNCNAME "CHM Elliptic Integral Function"

// I assume this shouldn't be here, but it's a static variable so hopefully
// it doesn't cause any problems...
static gsl_error_handler_t *old_error_handler=gsl_set_error_handler_off();

/* gsl_sf_expint_E1() isn't overridden by Boost so we'll
   compute its first ORDER derivatives by hand. The good news
   is that we only need to tell Boost that the first
   derivative is -exp(-x)/x. I didn't HAVE to do this to use this
   function in FunC (linear & cubic interpolation will work fine)
   but I wanted to demo the differentiation tables. */
// Define a local version of Boost's fvar
template <typename REAL_TYPE, size_t ORDER>
using CHM_fvar = boost::math::differentiation::autodiff_v1::detail::fvar<REAL_TYPE,ORDER>;

template <typename REAL_TYPE, size_t ORDER>
CHM_fvar<REAL_TYPE,ORDER> gsl_sf_expint_E1(CHM_fvar<REAL_TYPE,ORDER> const& cr){
  using boost::math::differentiation::make_fvar;
  // define the 0th derivative
  REAL_TYPE const d0 = gsl_sf_expint_E1(static_cast<REAL_TYPE>(cr));

  if(ORDER == 0)
    return d0;
  else{
    // Give Boost the first derivative of E1 and let it do the rest of the work
    auto x = make_fvar<REAL_TYPE, ORDER ? ORDER - 1 : 0>(static_cast<REAL_TYPE>(cr));
    auto const d1 = -exp(-x)/x; // d1 has the remaining ORDER-1 derivatives we need

    // define our function's taylor series expansion
    return cr.apply_coefficients(ORDER, [&d0, &d1](size_t i) { return i ? d1[i - 1] / i : d0; });
  }
}

template <typename T>
T MyFunction(T x){
  return -(log(x) + 0.577215 + gsl_sf_expint_E1(x));
}
