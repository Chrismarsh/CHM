/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "station.hpp"

#include "factory.hpp"

namespace pt = boost::property_tree;

typedef pt::ptree config_file;

class filter_base
{
public:
    filter_base( std::string name = "",
		 config_file input_cfg = pt::basic_ptree<std::string,std::string>() )
      : cfg(input_cfg),ID(name) {};
    virtual ~filter_base(){};

    virtual void init(std::shared_ptr<station>& station){};
    virtual void process(std::shared_ptr<station>& station){};

    bool is_nan(double variable)
    {
        if( variable == -9999.0)
            return true;

        if( std::isnan(variable) )
            return true;

        if( std::isinf(variable) )
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

/**
* Factory related convenience typedef and macros
*/
// Macros for easier registration of Filter implementations
// single argument ctor
typedef factory<filter_base, config_file> filter_factory;
#define REGISTER_FILTER_HPP(Implementation) \
private: \
   static const registration_helper<filter_base,Implementation, config_file> registrar;
#define STR_EXPAND(x) #x     // ensure x gets evaluated as a string,
#define STR(x) STR_EXPAND(x) // two-stage macro
#define REGISTER_FILTER_CPP(Implementation) \
   const registration_helper<filter_base,Implementation,config_file> Implementation::registrar(STR(Implementation));
