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
* \addtogroup modules
* @{
* \class PenmanMonteith_evaporation
* \brief Calculates Penman Monteith evaportation
*
* Calculates evapotranspiration via Penman-Monteith
*
* Depends:
* -  Solar shortwave "Qsi" [W/m^-1]
* -  Incoming Longwave radiation "Lin" [W/m^2]
* -  Albedo "albedo" [-]
* -  Saturated vapour pressure "es" [kpa]
* -  Actual vapour pressure "ea" [kpa]
*
* Provides:
* -
* - ET "ET" [mm/time]
*/
class PenmanMonteith_evaporation : public module_base
{
public:
    PenmanMonteith_evaporation();
    ~PenmanMonteith_evaporation();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

/**
@}
*/