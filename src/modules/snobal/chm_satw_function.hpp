#pragma once
#include "EvaluationFunctor.hpp"
#include <cmath>
#define SEA_LEVEL 1.013246e5
#define BOIL 3.7315e2
#define FUNC_W(x) (pow(1.e1, -7.90298 * (BOIL / x - 1.) + 5.02808 * log(BOIL / x) / log(1.e1) - 1.3816e-7 * (pow(1.e1, 1.1344e1 * (1. - x / BOIL)) - 1.) + 8.1328e-3 * (pow(1.e1, -3.49149 * (BOIL / x - 1.)) - 1.) + log(SEA_LEVEL) / log(1.e1)))
#define FUNCNAME_W "pow(1.e1, -7.90298 * (BOIL / x - 1.) + 5.02808 * log(BOIL / x) / log(1.e1) - 1.3816e-7 * (pow(1.e1, 1.1344e1 * (1. - x / BOIL)) - 1.) + 8.1328e-3 * (pow(1.e1, -3.49149 * (BOIL / x - 1.)) - 1.) + log(SEA_LEVEL) / log(1.e1))"

class MyFunction_satw final : public EvaluationFunctor<double,double>
{
public:
  double operator()(double x) override { return FUNC_W(x); }
  double deriv(double x) override
  {
  	return log(1.e1) * pow(1.e1, -7.90298 * (BOIL / x - 1.) + 5.02808 * log(BOIL / x) / log(1.e1) - 1.3816e-7 * (pow(1.e1, 1.1344e1 * (1. - x / BOIL)) - 1.) + 8.1328e-3 * (pow(1.e1, -3.49149 * (BOIL / x - 1.)) - 1.) + log(SEA_LEVEL) / log(1.e1)) * (7.90298 * BOIL / x / x - 5.02808 / log(1.e1) / x - 1.3816e-7 * log(1.e1)* pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x));
  }
  double deriv2(double x) override
  {
  	return log(1.e1) * pow(1.e1, -7.90298 * (BOIL / x - 1.) + 5.02808 * log(BOIL / x) / log(1.e1) - 1.3816e-7 * (pow(1.e1, 1.1344e1 * (1. - x / BOIL)) - 1.) + 8.1328e-3 * (pow(1.e1, -3.49149 * (BOIL / x - 1.)) - 1.) + log(SEA_LEVEL) / log(1.e1)) * (log(1.e1) * (7.90298 * BOIL / x / x - 5.02808 / log(1.e1) / x - 1.3816e-7 * log(1.e1)* pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x)) * (7.90298 * BOIL / x / x - 5.02808 / log(1.e1) / x - 1.3816e-7 * log(1.e1)* pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x)) + (-15.80596 * BOIL / x / x / x + 5.02808 / log(1.e1) / x / x - 1.3816e-7 * log(1.e1) * log(1.e1) * pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x) * (3.49149 * BOIL / x / x) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (-6.98298 * BOIL / x / x / x)));
  }
  double deriv3(double x) override
  {
  	return log(1.e1) * pow(1.e1, -7.90298 * (BOIL / x - 1.) + 5.02808 * log(BOIL / x) / log(1.e1) - 1.3816e-7 * (pow(1.e1, 1.1344e1 * (1. - x / BOIL)) - 1.) + 8.1328e-3 * (pow(1.e1, -3.49149 * (BOIL / x - 1.)) - 1.) + log(SEA_LEVEL) / log(1.e1)) * (log(1.e1) * log(1.e1) * (7.90298 * BOIL / x / x - 5.02808 / log(1.e1) / x - 1.3816e-7 * log(1.e1)* pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x)) * (7.90298 * BOIL / x / x - 5.02808 / log(1.e1) / x - 1.3816e-7 * log(1.e1)* pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x)) * (7.90298 * BOIL / x / x - 5.02808 / log(1.e1) / x - 1.3816e-7 * log(1.e1)* pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x)) + 3 * log(1.e1) * (7.90298 * BOIL / x / x - 5.02808 / log(1.e1) / x - 1.3816e-7 * log(1.e1)* pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x)) * (-15.80596 * BOIL / x / x / x + 5.02808 / log(1.e1) / x / x - 1.3816e-7 * log(1.e1) * log(1.e1) * pow(1.e1, 1.1344e1 * (1. - x / BOIL)) * (-1.1344e1 / BOIL) * (-1.1344e1 / BOIL) + 8.1328e-3 * log(1.e1) * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (3.49149 * BOIL / x / x) * (3.49149 * BOIL / x / x) + 8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (-6.98298 * BOIL / x / x / x)) + (8.1328e-3 * log(1.e1) * pow(1.e1, -3.49149 * (BOIL / x - 1.)) * (log(1.e1) * log(1.e1) *  (3.49149 * BOIL / x / x) * (3.49149 * BOIL / x / x) * (3.49149 * BOIL / x / x) + 3 * log(1.e1) * (3.49149 * BOIL / x / x) * (-6.98298 * BOIL / x / x / x) + (20.94894 * BOIL / x / x / x / x))));
  }
};
