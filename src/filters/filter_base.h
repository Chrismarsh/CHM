#pragma once

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "station.hpp"

namespace pt = boost::property_tree;

class filter_base
{
public:
    filter_base(){};
    virtual ~filter_base(){};
    virtual void process(boost::shared_ptr<station> station){};
    bool is_nan(double variable)
    {
        if( std::fabs(variable - -9999.0) < 1e-5)
            return true;
        if( std::isnan(variable) )
            return true;

        return false;
    }
    /**
     * Configuration file. If filter does not need one, then this will contain nothing
     */
    pt::ptree cfg;

    /**
    * ID of the module
    */
    std::string ID;
};