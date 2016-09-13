#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "interpolation.h"

#include <cstdlib>
#include <string>

#include <cmath>

#include <math.h>

/**
* \addtogroup modules
* @{
* \class const_llra_ta
* \brief Constant linear lapse rate adjustment.
*
* Constant linear lapse rate adjustment for air temperature of 0.0065 degC/m.
*
* Depends:
* - None
*
* Provides:
* - Air temperatue "t" [degC]
* - Air temperatue "const_llra_ta" [degC]
*/
class const_llra_ta : public module_base
{
public:
    const_llra_ta(config_file cfg);
    ~const_llra_ta();
    virtual void run(mesh_elem& face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };

};

/**
@}
*/