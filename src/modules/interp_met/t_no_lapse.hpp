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
* \class t_no_lapse
* \brief No adjustment
*
* Depends:
* - None
*
* Provides:
* - Air temperature "t" [degC]
* - Air temperatue "t_no_lapse" [degC]
*
*/
class t_no_lapse : public module_base
{
public:
    t_no_lapse(config_file cfg);
    ~t_no_lapse();
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