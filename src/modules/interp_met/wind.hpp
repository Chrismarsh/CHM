#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

/**
* \addtogroup modules
* @{
* \class Wind
* \brief Calculates Wind
*
* Calculates wind in a terrible way
*
* Depends:
* - Wind from met file "u" [m/s]
*
* Provides:
* - Wind "u" [m/s]
*/
class wind : public module_base
{
public:
    wind(std::string ID);
    ~wind();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

/**
@}
*/