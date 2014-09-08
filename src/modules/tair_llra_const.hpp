#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include "TPSpline.hpp"

#include <cstdlib>
#include <string>

#include <cmath>

#include <math.h>

class tair_llra_const : public module_base
{
public:
    tair_llra_const(std::string ID);
    ~tair_llra_const();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);



};
