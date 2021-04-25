#pragma once
#include <cmath>
#define SEA_LEVEL 1.013246e5
#define BOIL 3.7315e2
#define FUNC_W(x) (pow(1.e1, -7.90298 * (BOIL / x - 1.) + 5.02808 * log(BOIL / x) / log(1.e1) - 1.3816e-7 * (pow(1.e1, 1.1344e1 * (1. - x / BOIL)) - 1.) + 8.1328e-3 * (pow(1.e1, -3.49149 * (BOIL / x - 1.)) - 1.) + log(SEA_LEVEL) / log(1.e1)))
#define FUNCNAME_W "pow(1.e1, -7.90298 * (BOIL / x - 1.) + 5.02808 * log(BOIL / x) / log(1.e1) - 1.3816e-7 * (pow(1.e1, 1.1344e1 * (1. - x / BOIL)) - 1.) + 8.1328e-3 * (pow(1.e1, -3.49149 * (BOIL / x - 1.)) - 1.) + log(SEA_LEVEL) / log(1.e1))"

template <typename T>
T MyFunction_satw(T x){ return FUNC_W(x); }
