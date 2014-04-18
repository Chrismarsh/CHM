#pragma once

#include "../logger.h"
#include "module_base.hpp"
//#include "global.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

class Solar : public module_base
{
    public:
        Solar(std::string ID);
        ~Solar();
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);

};

