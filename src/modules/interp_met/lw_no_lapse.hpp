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
* \class lw_no_lapse
* \brief No longwave adjustment for elevation 9constant)
*
*
* Depends:
* - None
*
* Provides:
* - Longwave "ilwr" [W/m^2]
*
*/
class lw_no_lapse : public module_base
{
public:
    lw_no_lapse(config_file cfg);
    ~lw_no_lapse();
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
