#include "PhysConst.h"

namespace PhysConst
{
	const double Lv(const units::Celsius T)
	{
		return Lv() - 0.002361 * T.value;
	};
};
