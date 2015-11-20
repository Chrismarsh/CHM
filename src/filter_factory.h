#pragma once
#include "filter_base.h"
#include "logger.hpp"
#include "exception.hpp"
#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
namespace pt = boost::property_tree;

#include "macdonald_undercatch.h"
class filter_factory
{
public:
    filter_base* get(std::string ID, pt::ptree config);
};


