#pragma once

#include "module_base.hpp"
#include "logger.h"
#include "exception.hpp"

#include <string>


// Modules to build against
#include "TestModule.hpp"
#include "Solar.hpp"

class module_factory
{
public:
    module_base* get(std::string ID);
};
