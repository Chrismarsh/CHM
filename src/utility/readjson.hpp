#pragma once

#include <string>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

#include "jsonstrip.hpp"
#include "exception.hpp"

pt::ptree read_json(const std::string& path);