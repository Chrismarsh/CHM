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
* - Snow precip p_snow [m]
* - Liquid precip p_rain [m]
* - Fractional rain frac_precip_rain [-]
* - Fractional snow frac_precip_snow [-]
*/
class Harder_precip_phase : public module_base
{
public:
    Harder_precip_phase(config_file cfg);
    ~Harder_precip_phase();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);

    double b;
    double c;

};

/**
@}
*/