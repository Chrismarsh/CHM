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
* \class Longwave_from_obs
* \brief Annual longwave lapse rate adjustment to longwave from observations
*
* Annual longwave lapse rate adjustment to longwave from observations.
*
* Depends:
* - None
*
* Provides:
* - Longwave "ilwr" [W/m^2]
*
* Reference:
* Marty et al. (2002) (add full ref here)
*/
class Longwave_from_obs : public module_base
{
public:
    Longwave_from_obs(config_file cfg);
    ~Longwave_from_obs();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);
    virtual void init(mesh domain, boost::shared_ptr<global> global_param);
    struct data : public face_info
    {
        interpolation interp;
    };
};

/**
@}
*/
