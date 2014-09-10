#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

#include "triangulation.hpp"
#include "global.hpp"

/**
* \class module_base
* \brief Base class for individual modules
*
* A module must inherent from this class. This provides the interface for all module.
*/
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
    * Upon instantiation the module is assigned a value based upon its initialization order. Primairly used to dependency resolution. No real usage otherwise.
    */
    int IDnum;

    /**
    * Default constructor
    */
    module_base()
    {
        _provides = boost::make_shared<std::vector<std::string> >();
        _depends = boost::make_shared<std::vector<std::string> >();
        _depends_from_met = boost::make_shared<std::vector<std::string> >();
        IDnum = 0;
        ID = "uninitialized ID";
        _parallel_type = parallel::data; //default to data parallel, most common (?)
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
    * Needs to be implemented by each  data parallel module. This will be called and executed for each timestep
    * \param elem The terrain element (triangle) to be worked upon for an element parallel domain
    * \param global_param A pointer to the shared global paramter space with domain-wide paramters
    */
    virtual void run(mesh_elem &elem, boost::shared_ptr<global> global_param)
    {
    };

    /*
     * Needs to be implemented by each  domain parallel module. This will be called and executed for each timestep. Unique to domain parallel modules.
     * \param domain The terrain element (triangle) to be worked upon for an element parallel domain
     * \param global_parama A pointer to the shared global paramter space with domain-wide paramters
     */
    virtual void run(mesh domain, boost::shared_ptr<global> global_param)
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
     * List of the variables from other modules that this module depends upon
     */
    boost::shared_ptr<std::vector<std::string> > depends()
    {
        return _depends;
    }

    /**
    * List of the variables from the met files that this module depends upon
    */
    boost::shared_ptr<std::vector<std::string> > depends_from_met()
    {
        return _depends_from_met;
    }

protected:
    parallel _parallel_type;
    boost::shared_ptr<std::vector<std::string> > _provides;
    boost::shared_ptr<std::vector<std::string> > _depends;
    boost::shared_ptr<std::vector<std::string> > _depends_from_met;

};

/**
* Convenience typedef for modules.
*/
typedef boost::shared_ptr<module_base> module;