#pragma once
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <string>

#include "exception.hpp"
#include "logger.hpp"

namespace pt = boost::property_tree;

class config : public pt::ptree
{
public:
    config();
    ~config();
    void open(std::string file);

};


