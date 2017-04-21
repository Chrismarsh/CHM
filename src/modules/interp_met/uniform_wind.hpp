#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "math/coordinates.hpp"
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <cmath>


/**
* \addtogroup modules
* @{
* \class uniform_wind
* \brief Calculates wind speed and direction without any modification
*
* Depends:
* - Wind at reference height "U_R" [m/s]
* - Direction at reference height 'vw_dir' [degrees]
*
* Provides:
* - Wind "U_R" [m/s] at reference height
* - Wind direction 'vw_dir' [degrees]
*/
class uniform_wind : public module_base
{
public:
    uniform_wind(config_file cfg);
    ~uniform_wind();
    virtual void run(mesh domain);
    virtual void init(mesh domain);
    class lwinddata : public face_info
    {
    public:
        double curvature;
        interpolation interp;
        double corrected_theta;
        double W;
    };
};

/**
@}
*/
