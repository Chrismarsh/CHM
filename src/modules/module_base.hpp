#pragma once

#include <string>

class ModuleBase
{
public:
        std::string ID;
        virtual void run()=0;
        
        ModuleBase()
        {
            //nothing
        };
        virtual ~ModuleBase()
        {
            //nothing
        };
};

