//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
namespace pt = boost::property_tree;

#include "triangulation.hpp"
#include "global.hpp"
#include "timeseries/netcdf.hpp"
#include "factory.hpp"

//Create a process modules group in the doxygen docs to add individual modules to
/**
* \defgroup modules Process modules
*/

/**
* \class module_base
* \brief Base class for individual modules
*
* A module must inherent from this class. This provides the interface for all module.
*/

typedef pt::ptree config_file;

class module_base
{
public:

    /**
    * \enum parallel
    * A module must declare itself as either parallel data,
    * where the module only operates upon a single mesh element, or parallel domain, where
    * the module operates upon the domain in its entirety. The default is data parallel
    */
    enum class parallel
    {
        /**
        * Sets that this module is element parallel.
        * That is, this module can compute the solution at this element without needing to communicate with surrounding elements. There is no guarantee on element ordering
        */
                data,
        /**
        * Sets that this module is domain parallel.
        * That is, this module requires surrounding elements to compute its answer and that it is dependent upon the order of traversal of elements.
        */
                domain
    };
    /**
    * ID of the module
    */
    std::string ID;

    /**
     * Configuration file. If module does not need one, then this will contain nothing
     */
    pt::ptree cfg;

    /**
    * Upon instantiation the module is assigned a value based upon its initialization order. Primarily used to dependency resolution. No real usage otherwise.
    */
    int IDnum;

    /**
     * Global parameter store
     */
    boost::shared_ptr<global> global_param;

    /**
    * Default constructor
    */
    module_base(){};

    /**
     * Consturctor that initializes everything
     */
    module_base(std::string name = "",
		parallel type = parallel::data,
		config_file input_cfg = pt::basic_ptree<std::string,std::string>())
      :    ID(name), cfg(input_cfg), IDnum(0),_parallel_type(type)
    {
        _provides = boost::make_shared<std::vector<std::string> >();
        _depends = boost::make_shared<std::vector<std::string> >();
        _depends_from_met = boost::make_shared<std::vector<std::string> >();
        _optional = boost::make_shared<std::vector<std::string> >();
        _conflicts = boost::make_shared<std::vector<std::string> >(); //modules that we explicitly cannot be run alongside. Use sparingly
        global_param = nullptr;

        //nothing
    };

    /**
    * Default destructor
    */
    virtual ~module_base()
    {
        //nothing
    };

    /**
     * Checkpoint (save state) the current module. By default this errors out if the modules does not support checkpointing.
     * @param domain
     * @param data
     */
    virtual void checkpoint(mesh& domain, netcdf& chkpt)
    {

      //TODO: Add default check for the assumption that module does not support serialization
    };

    virtual void load_checkpoint(mesh& domain, netcdf& chkpt)
    {

        //TODO: Add default check for the assumption that module does not support serialization
    };

    /**
    * Needs to be implemented by each  data parallel module. This will be called and executed for each timestep
    * \param face The terrain element (triangle) to be worked upon for an element parallel domain
    * \param global_param A pointer to the shared global paramter space with domain-wide paramters
    */
    virtual void run(elem&&face)
    {
    };

    /*
     * Needs to be implemented by each  domain parallel module. This will be called and executed for each timestep. Unique to domain parallel modules.
     * \param domain The entier terrain mesh
     * \param global_parama A pointer to the shared global paramter space with domain-wide paramters
     */
    virtual void run(mesh& domain)
    {
    };

    /*
     * Optional function to run after the dependency constructor call, but before the run function is called. Used to perform any initalization.
     * \param domain The entire terrain mesh
     */
    virtual void init(mesh& domain)
    {

    };

    /*
     * Returns the module's parallel type
     * \return the parallel type
     */
    parallel parallel_type()
    {
        return _parallel_type;
    }

    /**
    * List of the variables that this module provides.
    */
    boost::shared_ptr<std::vector<std::string> > provides()
    {
        return _provides;
    }

    /**
     * Set a variable that this module provides
     */
    void provides(const std::string& variable)
    {
        if(variable.find_first_of("\t ") != std::string::npos)
            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Variable " + variable +" has a space. This is not allowed."));

        _provides->push_back(variable);
    }

    /**
     * List of the variables from other modules that this module depends upon
     */
    boost::shared_ptr<std::vector<std::string> > depends()
    {
        return _depends;
    }

    /**
    * Modules we conflict with and absolutely cannot run alongside. Use sparingly.
    */
    void conflicts(const std::string& variable)
    {
        if(variable.find_first_of("\t ") != std::string::npos)
            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Variable " + variable +" has a space. This is not allowed."));

        _conflicts->push_back(variable);
    }

    boost::shared_ptr<std::vector<std::string> > conflicts()
    {
        return _conflicts;
    }


    /**
     * List of the optional depends variables from other modules that this module depends upon
     */
    boost::shared_ptr<std::vector<std::string> > optionals()
    {
        return _optional;
    }


    /**
     * Set a variable, from another module, that this module depends upon
     */
    void depends(const std::string& variable)
    {
        _depends->push_back(variable);
    }

    /**
    * List of the variables from the met files that this module depends upon
    */
    boost::shared_ptr<std::vector<std::string> > depends_from_met()
    {
        return _depends_from_met;
    }

    /**
     * Set a variable, from a met file, that this module depends upon
     */
    void depends_from_met(const std::string& variable)
    {
        _depends_from_met->push_back(variable);
    }

    /**
     * Set an optional (not required) variable, from another module, that this module depends upon.
     *
     */
    void optional(const std::string& variable)
    {
        _optional->push_back(variable);
        _optional_found.insert( std::pair<std::string,bool>(variable,false));
    }

    /**
     * Checks if an optional variable was found
     */
    bool has_optional(const std::string& variable)
    {
        auto it = _optional_found.find(variable);

        //asked for a variable that isn't optional, just return false
        if(it ==_optional_found.end())
            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Requested a non-optional variable"));
//            return false;

        return it->second; //return if we found it


    }

    /**
     * If you want to skip evaluating this current face, call this to set all provides outputs to nan
     * E.g., called if the is_water, is_glacier, etc is true
     */
    void set_all_nan_on_skip(elem& face)
    {
        for(auto& itr: *_provides)
        {
            face[itr]=-9999.;
        }
    }
    /**
     * Set that an optional variable was found
     */
    void set_optional_found(const std::string& variable)
    {
        _optional_found[variable]=true;
    }

    bool is_nan(const double& variable)
    {
        if( std::fabs(variable - -9999.0) < 1e-5)
            return true;
        if( std::isnan(variable) )
            return true;

        return false;
    }

    bool is_water(elem& face)
    {
        bool is = false;

        if(face.has_parameter("landcover"))
        {
            int LC = face.get_parameter("landcover");
            is = global_param->parameters.get<bool>("landcover." + std::to_string(LC) + ".is_water",false);
        }
        return is;
    }

    bool is_glacier(elem& face)
    {
        bool is = false;

        if(face.has_parameter("landcover"))
        {
            int LC = face.get_parameter("landcover");
            is = global_param->parameters.get<bool>("landcover." + std::to_string(LC) + ".is_glacier",false);
        }
        return is;
    }

protected:
    parallel _parallel_type;
    boost::shared_ptr<std::vector<std::string> > _provides;
    boost::shared_ptr<std::vector<std::string> > _depends;
    boost::shared_ptr<std::vector<std::string> > _depends_from_met;
    boost::shared_ptr<std::vector<std::string> > _optional;
    boost::shared_ptr<std::vector<std::string> > _conflicts;



    //lists the options that were found
    std::map<std::string,bool> _optional_found;


};

/**
* Convenience typedef for modules.
*/
typedef boost::shared_ptr<module_base> module;

/**
* Factory related convenience typedef and macros
*/
// Macros for easier registration of Module implementations
// single argument ctor
typedef factory<module_base, config_file> module_factory;
#define REGISTER_MODULE_HPP(Implementation) \
private: \
   static const registration_helper<module_base,Implementation,config_file> registrar;
#define STR_EXPAND(x) #x     // ensure x gets evaluated as a string,
#define STR(x) STR_EXPAND(x) // two-stage macro
#define REGISTER_MODULE_CPP(Implementation) \
   const registration_helper<module_base,Implementation,config_file> Implementation::registrar(STR(Implementation));
