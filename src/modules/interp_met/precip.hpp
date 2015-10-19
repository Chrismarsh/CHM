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
* \class Precip
* \brief Calculates precip
*
* Spatially distributes precip via spline
*
* Depends:
* - Precip from met file "p" [mm]
*
* Provides:
* - Precip "p" [m/s]
*/
class precip : public module_base
{
public:
    precip();
    ~precip();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

/**
@}
*/