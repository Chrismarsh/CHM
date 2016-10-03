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

//#include <meteoio/MeteoIO.h>
/**
* \addtogroup modules
* @{
* \class Harder_precip_phase
* \brief Calculates precip phase based
*
* Calculates precipitation phase via falling hydrometeor energy balance
*
* Depends:
* - Air temperature "t" [C]
* - Relative Humidity 'rh' [C]
* - Precip "p" [mm]
*
* Provides:
* - Snow precip p_snow [mm]
* - Liquid precip p_rain [mm]
* - Fractional rain frac_precip_rain [-]
* - Fractional snow frac_precip_snow [-]
*/
class Harder_precip_phase : public module_base
{
public:
    Harder_precip_phase(config_file cfg);
    ~Harder_precip_phase();
    virtual void run(mesh_elem& face);

    double b;
    double c;

};

/**
@}
*/