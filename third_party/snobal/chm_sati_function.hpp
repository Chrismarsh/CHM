#pragma once
#include <cmath>
#define FREEZE 2.7316e2
#define FUNC_I(x) (pow(1.e1, -9.09718 * ((FREEZE / x) - 1.) - 3.56654 * log(FREEZE / x) / log(1.e1) + 8.76793e-1 * (1. - (x / FREEZE)) + log(6.1071) / log(1.e1)) * 100.0)
#define FUNCNAME_I "pow(1.e1, -9.09718 * ((FREEZE / x) - 1.) - 3.56654 * log(FREEZE / x) / log(1.e1) + 8.76793e-1 * (1. - (x / FREEZE)) + log(6.1071) / log(1.e1)) * 100.0"

template <typename T>
T MyFunction_sati(T x){ return FUNC_I(x); }
