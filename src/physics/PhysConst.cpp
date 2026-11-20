#include "PhysConst.h"

namespace PhysConst
{
        double Lv(const double T) {
            return Lv() - 0.002361 * T;
	};
};
