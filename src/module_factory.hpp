#pragma once

#include "module_base.hpp"
#include "logger.h"
#include <string>


// Modules to build against
#include "TestModule.hpp"

class ModuleFactory
{
public:
    ModuleBase* get(std::string ID);
};
