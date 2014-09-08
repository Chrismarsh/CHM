#pragma once


#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include "TPSpline.hpp"

#include <cstdlib>
#include <string>

#include <cmath>

#include <math.h>

class rh_llra_var : public module_base
{
public:
    rh_llra_var(std::string ID);
    ~rh_llra_var();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};
