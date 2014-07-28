#pragma once

#include "module_base.hpp"
#include "logger.h"
#include "exception.hpp"

#include <string>


// Modules to build against
#include "solar.hpp"
#include "terrain_shadows.hpp"
class module_factory
{
public:
    module_base* get(std::string ID);
};
