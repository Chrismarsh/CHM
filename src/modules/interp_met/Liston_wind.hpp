#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "math/distance.hpp"
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <cmath>


/**
* \addtogroup modules
* @{
* \class Liston_wind
* \brief Calculates wind speed and direction following Liston 2006
*
* Calculates windspeeds using a terrain curvature following Liston 2006.
* Defaults to 300 m distance. Configurable using "distance".
 * Liston, G. E., & Elder, K. (2006). A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). Journal of Hydrometeorology, 7(2), 217–234. http://doi.org/10.1175/JHM486.1
* Depends:
* - Wind at reference height "U_R" [m/s]
* - Direction at reference height 'vw_dir' [degrees]
*
* Provides:
* - Wind "U_R" [m/s] at reference height
* - Wind direction 'vw_dir' [degrees]
*/
class Liston_wind : public module_base
{
public:
    Liston_wind(config_file cfg);
    ~Liston_wind();
    virtual void run(mesh domain, boost::shared_ptr<global> global_param);
    virtual void init(mesh domain, boost::shared_ptr<global> global_param);
    class lwinddata : public face_info
    {
    public:
        double curvature;
        interpolation interp;
        double corrected_theta;
        double W;
    };
    double distance;
};

/**
@}
*/