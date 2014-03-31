#pragma once

#include "../logger.h"
#include "module_base.hpp"
#include <string>

class TestModule : public ModuleBase
{
    public:
        TestModule(std::string ID);
        ~TestModule();
        void run();

};

