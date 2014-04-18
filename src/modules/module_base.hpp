#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

#include "triangle.h"
#include "global.hpp"


typedef triangle mesh_elem;

class module_base
{
public:
        std::string ID;
        
        module_base()
        {
            //nothing
        };
        virtual ~module_base()
        {
            //nothing
        };
        
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param)=0;
};



typedef boost::shared_ptr< module_base > module;