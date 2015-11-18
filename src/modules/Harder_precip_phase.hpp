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
* \brief Calculates Thornton_p phase based
*
* Calculates wind in a terrible way
*
* Depends:
* - Air temperature "t" [degrees]
* - Relative Humidity 'rh' [degrees]
*
* Provides:
* - Snow Thornton_p p_snow [m]
* - Liquid Thornton_p p_rain [m]
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