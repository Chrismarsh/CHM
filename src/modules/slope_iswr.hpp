#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <meteoio/MeteoIO.h>
#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

/**
* \addtogroup modules
* @{
* \class slope_iswr
* \brief Shortwave radiation on a slope
*
* Calculates incoming direct-beam shortwave solar radiation to slope
*
* Depends:
* - Shortwave all beam "iswr" [W/m^2]
* - Shortwave direct "iswr_direct" [W/m^2]
* - Shortwave diffuse "iswr_diffuse" [W/m^2]

* Provides:
* - Shortwave all beam "iswr" [W/m^2]
* - Shortwave direct "iswr_direct" [W/m^2]
* - Shortwave diffuse "iswr_diffuse" [W/m^2]
*/
class slope_iswr : public module_base
{
    public:
        slope_iswr(config_file cfg);
        ~slope_iswr();
        virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


        bool assume_no_slope;
};

/**
@}
*/