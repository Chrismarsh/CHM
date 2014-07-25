#pragma once

#include "../logger.h"
#include "triangulation.h"
#include "module_base.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

class terrain_shadow : public module_base
{
    public:
        terrain_shadow(std::string ID);
        ~terrain_shadow();
        virtual void run(mesh domain, boost::shared_ptr<global> global_param);
        

};

