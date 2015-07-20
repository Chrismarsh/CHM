#pragma once

#include "module_base.hpp"
#include "logger.hpp"
#include "exception.hpp"

#include <string>


// Modules to build against
#include "solar.hpp"
#include "Marsh_shading_iswr.hpp"
#include "const_llra_ta.hpp"
#include "Liston_monthly_llra_rh.hpp"
#include "Liston_monthly_llra_ta.hpp"
#include "Sicart_ilwr.hpp"
#include "wind.hpp"
#include "PenmanMonteith_evaporation.hpp"
#include "precip.hpp"


class module_factory
{
public:
    module_base* get(std::string ID);
};
