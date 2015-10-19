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
* \class Sicart_ilwr
* \brief Calculates longwave radiation from the atmosphere (no skyview correction) for clear and cloudy days.
*
* Calculates evapotranspiration via Penman-Monteith
*
* Depends:
* -  Transmittance "atm_trans" [-]
* -  Actual vapour pressure "ea" [kpa]
* -  Daily atmospheric transmittance "atm_trans" [-]
*
* Provides:
* - Incoming longwave  "Lin" [W/m^2]
*/
class Sicart_ilwr : public module_base
{
public:
    Sicart_ilwr();
    ~Sicart_ilwr();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};

/**
@}
*/