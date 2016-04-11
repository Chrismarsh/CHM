#pragma once

#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

#include "jsonstrip.cpp"

#include "exception.hpp"

pt::ptree read_json(std::string path);