#pragma once

#include "../logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <cstdlib>
#include <string>

#include <cmath>

#include <math.h>

/**
* \addtogroup modules
* @{
* \class t_monthly_lapse
* \brief Constant monthly linear lapse rate adjustment
*
* Constant monthly linear lapse rate adjustment for air temperature.
* Taken as input from user or default values for Canadian Rockies from:
 * *
* Joseph M. Shea, Shawn J. Marshall and Joanne M. Livingston
* Arctic, Antarctic, and Alpine Research
* Vol. 36, No. 2 (May, 2004), pp. 272-279
*
* Depends:
* - None
*
* Provides:
* - Air temperature "t" [degC]
* - Air temperatue "t_monthly_lapse" [degC]
*

*/
class t_monthly_lapse : public module_base
{
public:
    t_monthly_lapse(config_file cfg);
    ~t_monthly_lapse();
    virtual void run(mesh_elem& face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };
    double MLR[12];
};

/**
@}
*/