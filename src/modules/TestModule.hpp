#pragma once

#include "../logger.h"
#include "module_base.hpp"

#include <cstdlib>
#include <string>

class TestModule : public module_base
{
    public:
        TestModule(std::string ID);
        ~TestModule();
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);

};

