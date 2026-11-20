#pragma once

#include <functional>
#include <cmath>
#include <limits>
namespace math
{
    enum class Root { Found, DidNotConverge, NoRootGuaranteed };
    struct Result
    {
        Root root;
        double value;
        size_t iterations;
    };

    inline Result bisection(
        std::function<double(double)> f,
        double low_guess,
        double  high_guess,
        double tol = 1e-10,
        size_t max_iter = 1000
    ) {
        double fa = f(low_guess), fb = f(high_guess);

        if (fa * fb > 0) 
            return Result {
                Root::NoRootGuaranteed,
                std::numeric_limits<double>::quiet_NaN(),
                0};

        for (size_t i = 0; i < max_iter; ++i) {
            double c = low_guess + (high_guess - low_guess) / 2;
            double fc = f(c);

            if (std::abs(fc) < tol || (high_guess - low_guess) / 2 < tol) 
                return Result{Root::Found,c,i} ;

            if (fa * fc < 0) {
                high_guess = c;
                fb = fc;
            } else {
                low_guess = c;
                fa = fc;
            }
        }

        return Result{Root::DidNotConverge,
            std::numeric_limits<double>::quiet_NaN(),
            max_iter};  // Did not converge
    };
};
