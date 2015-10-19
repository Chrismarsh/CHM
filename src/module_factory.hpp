#pragma once

#include "module_base.hpp"
#include "logger.hpp"
#include "exception.hpp"
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
namespace pt = boost::property_tree;

// Modules to build against
#include "solar.hpp"
#include "Marsh_shading_iswr.hpp"
#include "const_llra_ta.hpp"
#include "Liston_monthly_llra_rh.hpp"
#include "Liston_monthly_llra_ta.hpp"
#include "Sicart_ilwr.hpp"
#include "Liston_wind.hpp"
#include "PenmanMonteith_evaporation.hpp"
#include "precip.hpp"
#include "Walcek_atm_trans.hpp"
#include "Harder_precip_phase.hpp"

namespace pt = boost::property_tree;
class module_factory
{
public:
    module_base* get(std::string ID, pt::ptree& cfg);
};
