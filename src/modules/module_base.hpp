#pragma once

#include <string>
#include <boost/shared_ptr.hpp>

#include "triangle.h"

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
        
        virtual void run(mesh_elem& elem)=0;
};



typedef boost::shared_ptr< module_base > module;