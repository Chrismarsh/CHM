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
#include <boost/math/tools/roots.hpp>

/**
* \addtogroup modules
* @{
* \class Harder_precip_phase
* \brief Calculates precip phase based
*
* Calculates wind in a terrible way
*
* Depends:
* - Air temperature "t" [degrees]
* - Relative Humidity 'rh' [degrees]
*
* Provides:
* - Snow precip p_snow [m]
* - Liquid precip p_rain [m]
*/
class Harder_precip_phase : public module_base
{
public:
    Harder_precip_phase();
    ~Harder_precip_phase();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

/**
@}
*/