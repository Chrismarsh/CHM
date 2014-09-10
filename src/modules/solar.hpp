#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

/**
* \class Solar
* \breif Calculates shortwave radiation
*
* Calculates incoming direct-beam shortwave solar radiation, with no correction for transmissivity.
*/
class Solar : public module_base
{
    public:
        Solar(std::string ID);
        ~Solar();
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

