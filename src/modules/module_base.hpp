#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

#include "triangulation.h"
#include "global.hpp"

/*
 * The base class for individual modules
 */
class module_base
{
public:
    enum class parallel
    {
        data,
        domain
    };
        //Identification string for this module
        std::string ID;
        int IDnum;
        
        /*
         * Empty constructor
         */
        module_base()
        {
            _provides = boost::make_shared<std::vector<std::string> >();
            _depends = boost::make_shared<std::vector<std::string> >();
            IDnum = 0;
            ID = "uninitialized";
            _parallel_type = parallel::data; //default to data parallel, most common (?)
            //nothing
        };
        virtual ~module_base()
        {
            //nothing
        };
        
        /*
         * Needs to be implimented by each module. This will be called and executed for each timestep
         * @param elem The terrain element (triangle) to be worked upon for an element parallel domain
         * @param global_parama A pointer to the shared global paramter space with domain-wide paramters
         */
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param){};
        
        /*
         * Needs to be implimented by each module. This will be called and executed for each timestep. Unique to domain parallel modules;
         * @param elem The terrain element (triangle) to be worked upon for an element parallel domain
         * @param global_parama A pointer to the shared global paramter space with domain-wide paramters
         */
        virtual void run(mesh domain, boost::shared_ptr<global> global_param){};
        
        /*
         * Sets that this module is element parallel. That is, this module can compute the solution at this element without needing to communicate with surrounding elements.
         * Elem_parallel modules will be run over the domain in an unspecified order. There is no guarantee on element ordering
         */
      
        /*
         * Sets that this module is domain parallel. That is, this module requires surrounding elements to compute it's answer and that it is dependent upon the order of traversal of elements.
         */

        /*
         * Returns if module is domain parallel
         * @return true if module is domain parallel 
         */
        parallel parallel_type()
        {
            return _parallel_type;
        }

        
        boost::shared_ptr<std::vector<std::string> > provides()
        {
            return _provides;
        }
        boost::shared_ptr<std::vector<std::string> > depends()
        {
            return _depends;
        }
protected:
        parallel _parallel_type;
        boost::shared_ptr<std::vector<std::string> > _provides;
        boost::shared_ptr<std::vector<std::string> > _depends;
        
};

typedef boost::shared_ptr< module_base > module;