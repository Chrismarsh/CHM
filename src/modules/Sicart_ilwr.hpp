#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#include "math/coordinates.hpp"
#define _USE_MATH_DEFINES
#include <math.h>
#include <meteoio/MeteoIO.h>
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
* -  Daily atmospheric transmittance "atm_trans" [-]
*
* Provides:
* - Incoming longwave  "ilwr" [W/m^2]
*/
class Sicart_ilwr : public module_base
{
public:
    Sicart_ilwr(config_file cfg);
    ~Sicart_ilwr();
    virtual void run(mesh_elem& face);
    void init(mesh domain);


};

/**
@}
*/