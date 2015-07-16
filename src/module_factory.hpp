#pragma once

#include "module_base.hpp"
#include "logger.hpp"
#include "exception.hpp"

#include <string>


// Modules to build against
#include "solar.hpp"
#include "terrain_shadows.hpp"
#include "tair_llra_const.hpp"
#include "rh_llra_var.hpp"
#include "tair_llra_lookup.hpp"
#include "atm_trans_annandale.hpp"
#include "longwave_sicart.hpp"
#include "wind.hpp"
#include "evap_penman_monteith.h"
#include "precip.hpp"
#include "leaky_bucket.hpp"

class module_factory
{
public:
    module_base* get(std::string ID);
};
