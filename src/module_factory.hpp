#pragma once

#include "module_base.hpp"
#include "logger.h"
#include "exception.hpp"

#include <string>


// Modules to build against
#include "TestModule.hpp"

class module_factory
{
public:
    ModuleBase* get(std::string ID);
};
