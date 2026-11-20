#include "PhysConst.hpp"
#include "PhysFunc.hpp"

namespace PhysFunc
{
    double Lv(const double T) {
        return PhysConst::Lv - 0.002361 * T;
    };
};
