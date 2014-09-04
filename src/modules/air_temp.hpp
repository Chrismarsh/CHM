#pragma once

#include "../logger.h"
#include "triangulation.h"
#include "module_base.hpp"

#include "TPSpline.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

class air_temp : public module_base
{
public:
    air_temp(std::string ID);
    ~air_temp();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);



};
