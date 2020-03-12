#pragma once
#include <func/EvaluationFunctor.hpp>
#include <cmath>
#define FREEZE 2.7316e2
#define FUNC_I(x) (pow(1.e1, -9.09718 * ((FREEZE / x) - 1.) - 3.56654 * log(FREEZE / x) / log(1.e1) + 8.76793e-1 * (1. - (x / FREEZE)) + log(6.1071) / log(1.e1)) * 100.0)
#define FUNCNAME_I "pow(1.e1, -9.09718 * ((FREEZE / x) - 1.) - 3.56654 * log(FREEZE / x) / log(1.e1) + 8.76793e-1 * (1. - (x / FREEZE)) + log(6.1071) / log(1.e1)) * 100.0"

class MyFunction_sati final : public EvaluationFunctor<double,double>
{
public:
  double operator()(double x) override { return FUNC_I(x); }
  double deriv(double x) override
  {
  	return 100.0 * log(1.e1) * pow(1.e1, -9.09718 * ((FREEZE / x) - 1.) - 3.56654 * log(FREEZE / x) / log(1.e1) + 8.76793e-1 * (1. - (x / FREEZE)) + log(6.1071) / log(1.e1)) * (9.09718 * FREEZE / x / x + 3.56654 / x / log(1.e1) - 8.76793e-1 / FREEZE);
  }
  double deriv2(double x) override
  {
  	return 100.0 * log(1.e1) * pow(1.e1, -9.09718 * ((FREEZE / x) - 1.) - 3.56654 * log(FREEZE / x) / log(1.e1) + 8.76793e-1 * (1. - (x / FREEZE)) + log(6.1071) / log(1.e1)) * ((9.09718 * FREEZE / x / x + 3.56654 / x / log(1.e1) - 8.76793e-1 / FREEZE) * (9.09718 * FREEZE / x / x + 3.56654 / x / log(1.e1) - 8.76793e-1 / FREEZE) * log(1.e1) + (-18.19436 * FREEZE / x / x / x - 3.56654 / x / x / log(1.e1)));
  }
  double deriv3(double x) override
  {
  	return 100.0 * log(1.e1) * pow(1.e1, -9.09718 * ((FREEZE / x) - 1.) - 3.56654 * log(FREEZE / x) / log(1.e1) + 8.76793e-1 * (1. - (x / FREEZE)) + log(6.1071) / log(1.e1)) * ((9.09718 * FREEZE / x / x + 3.56654 / x / log(1.e1) - 8.76793e-1 / FREEZE) * (9.09718 * FREEZE / x / x + 3.56654 / x / log(1.e1) - 8.76793e-1 / FREEZE) * (9.09718 * FREEZE / x / x + 3.56654 / x / log(1.e1) - 8.76793e-1 / FREEZE) * log(1.e1) * log(1.e1) + 3.0 * log(1.e1) * (9.09718 * FREEZE / x / x + 3.56654 / x / log(1.e1) - 8.76793e-1 / FREEZE) * (-18.19436 * FREEZE / x / x / x - 3.56654 / x / x / log(1.e1)) + (54.58308 * FREEZE / x / x / x / x + 7.13308 / x / x / x / log(1.e1)));
  }
};
