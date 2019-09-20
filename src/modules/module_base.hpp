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

/**
     * \enum SpatialType
     * Module "depends" variables must specify where the spatial extent the look at for their dependent variables.
     *  local     = only this face element
     *  neighbour = this and nearest-neighbour face elements
     *  distance  = this and all face elements within a fixed distance
     */
    enum class SpatialType
    {

        /**
         * Sets that dependent variable is local.
         * That is, this module requires the  compute the solution at this element without needing to communicate with
         * surrounding elements. There is no guarantee on element ordering
         */
        local,
        /**
         * Sets that this module is domain parallel.
         * That is, this module requires surrounding elements to compute its answer and that it is dependent upon the
         * order of traversal of elements.
         */
        neighbour,
        /**
         * Sets that this module is domain parallel.
         * That is, this module requires surrounding elements to compute its answer and that it is dependent upon the
         * order of traversal of elements.
         */
        distance
    };
    /**
     * Comprehensive info for spatial variable dependencies

     */
    struct variable_info
    {
        bool is_distance_set;
        SpatialType spatial_type;
        double spatial_distance;
        std::string name;
      // no need for default constructor
        variable_info() = delete;
      // constructing by name only assumes local
        variable_info(std::string name) : spatial_type{SpatialType::local}, name{name} {}
      // constructing by spatial type explicitly
      //  - if type is distance, need to make sure distance gets set before usage (useful for setting via options)
        variable_info(std::string name, SpatialType st) : spatial_type{st}, name{name}
        {
            switch (st)
            {
            case SpatialType::distance:
                is_distance_set = false;
            }
        }
      // explicit construction with distance intended for modules with a fixed distance
        variable_info(std::string name, SpatialType st, double distance) : name{name}, spatial_type{st}
        {
            switch (st)
            {
            case SpatialType::distance:
                is_distance_set = true;
                spatial_distance = distance;
            default:
                BOOST_THROW_EXCEPTION(
                    module_error() << errstr_info(
                        "Distance specified with local or neighbour SpatialType. This is not allowed."));
            }
        }
    };


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
        _provides = boost::make_shared<std::vector<variable_info>>();
        _provides_parameters = boost::make_shared<std::vector<std::string>>();
        _depends = boost::make_shared<std::vector<variable_info>>();
        _depends_from_met = boost::make_shared<std::vector<std::string>>();
        _optional = boost::make_shared<std::vector<std::string>>();
        _conflicts = boost::make_shared<std::vector<std::string>>(); // modules that we explicitly cannot be run
                                                                     // alongside. Use sparingly
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
    virtual void run(mesh_elem&face)
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
    boost::shared_ptr<std::vector<variable_info> > provides()
    {
        return _provides;
    }

    /**
    * List of the parameters that this module provides.
    */
    boost::shared_ptr<std::vector<std::string> > provides_parameter()
    {
        return _provides_parameters;
    }

    /**
     * Set a variable that this module provides
     */
    void provides(const std::string& name)
    {
        if(name.find_first_of("\t ") != std::string::npos)
            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Variable " + name +" has a space. This is not allowed."));

        _provides->push_back(variable_info(name));
    }

    /**
     * Set a variable that this module provides
     */
    void provides(const std::string& name, SpatialType st)
    {
        if(name.find_first_of("\t ") != std::string::npos)
            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Variable " + name +" has a space. This is not allowed."));

        _provides->push_back(variable_info(name, st));
    }

    /**
     * Set a variable that this module provides
     */
    void provides(const std::string& name, SpatialType st, double distance)
    {
        if(name.find_first_of("\t ") != std::string::npos)
            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Variable " + name +" has a space. This is not allowed."));

        _provides->push_back(variable_info(name, st, distance));
    }


    /**
     * Set a parameter that this module provides
     */
    void provides_parameter(const std::string& variable)
    {
        if(variable.find_first_of("\t ") != std::string::npos)
            BOOST_THROW_EXCEPTION(module_error() << errstr_info ("Variable " + variable +" has a space. This is not allowed."));

        _provides_parameters->push_back(variable);
    }

    /**
     * List of the variables from other modules that this module depends upon
     */
    boost::shared_ptr<std::vector<variable_info>> depends() { return _depends; }

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
    void depends(const std::string& name) { _depends->push_back(variable_info(name)); }
    void depends(const std::string& name, SpatialType st) { _depends->push_back(variable_info(name,st)); }
    void depends(const std::string& name, SpatialType st, double distance) { _depends->push_back(variable_info(name,st,distance)); }

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
    void set_all_nan_on_skip(mesh_elem& face)
    {
        for(auto& itr: *_provides)
        {
            (*face)[itr.name]=-9999.;
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

    bool is_water(mesh_elem& face)
    {
        bool is = false;

        if(face->has_parameter("landcover"_s))
        {
            int LC = face->parameter("landcover"_s);
            is = global_param->parameters.get<bool>("landcover." + std::to_string(LC) + ".is_water",false);
        }
        return is;
    }

    bool is_glacier(mesh_elem& face)
    {
        bool is = false;

        if(face->has_parameter("landcover"_s))
        {
            int LC = face->parameter("landcover"_s);
            is = global_param->parameters.get<bool>("landcover." + std::to_string(LC) + ".is_glacier",false);
        }
        return is;
    }

protected:
    parallel _parallel_type;
    boost::shared_ptr<std::vector<variable_info>> _provides;
    boost::shared_ptr<std::vector<std::string>> _provides_parameters;
    boost::shared_ptr<std::vector<variable_info>> _depends;
    boost::shared_ptr<std::vector<std::string>> _depends_from_met;
    boost::shared_ptr<std::vector<std::string>> _optional;
    boost::shared_ptr<std::vector<std::string>> _conflicts;

    // lists the options that were found
    std::map<std::string, bool> _optional_found;
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
