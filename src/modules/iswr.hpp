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
* \class iswr
* \brief
*
* This aggregates the _no_slope direct and diffuse flat plane estimates provided by other modules and calculates a unified
 * variant that includes (if required) slope effects and applies shadow masks from shadowing modules.
*
* Depends:
* - Shortwave direct on a flat plane "iswr_direct_no_slope" [W/m^2]
* - Shortwave diffuse on a flat plane "iswr_diffuse_no_slope" [W/m^2]

* Provides:
* - Shortwave all beam "iswr" [W/m^2]
* - Shortwave direct "iswr_direct" [W/m^2]
* - Shortwave diffuse "iswr_diffuse" [W/m^2]
*/
class iswr : public module_base
{
    public:
        iswr(config_file cfg);
        ~iswr();
        virtual void run(mesh_elem& face);


        bool assume_no_slope;

        // if we are using obs, then our obs implicitily have a cosine correction.
        // This needs to be undo prior to the correction for slope
        bool already_cosine_corrected;

};

/**
@}
*/