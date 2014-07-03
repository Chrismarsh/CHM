#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

#include "triangle.h"
#include "global.hpp"


typedef triangle mesh_elem;

/*
 * The base class for individual modules
 */
class module_base
{
public:
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
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param)=0;
        
        /*
         * Sets that this module is element parallel. That is, this module can compute the solution at this element without needing to communicate with surrounding elements.
         * Elem_parallel modules will be run over the domain in an unspecified order. There is no guarantee on element ordering
         */
        void set_is_elem_parallel()
        {
            _is_elem_parallel = true;
            _is_domain_parallel = false;
        }
        
        /*
         * Sets that this module is domain parallel. That is, this module requires surrounding elements to compute it's answer and that it is dependent upon the order of traversal of elements.
         */
        void set_is_domain_parallel()
        {
            _is_elem_parallel = false;
            _is_domain_parallel = true;
        }
        
        /*
         * Returns if module is domain parallel
         * @return true if module is domain parallel 
         */
        bool is_domain_parallel()
        {
            return _is_domain_parallel;
        }
        
        /*
        * Returns if module is element parallel
        * @return true if module is element parallel 
        */
        bool is_elem_parallel()
        {
            return _is_elem_parallel;
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
        bool _is_elem_parallel;
        bool _is_domain_parallel;
        boost::shared_ptr<std::vector<std::string> > _provides;
        boost::shared_ptr<std::vector<std::string> > _depends;
        
};

typedef boost::shared_ptr< module_base > module;