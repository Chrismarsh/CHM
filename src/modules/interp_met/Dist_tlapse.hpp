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
* \class Dist_tlapse
* \brief Distributed lapse rate from forcing file, changes per timestep
*
* 
*
* Depends:
* - Air temperature "t" [degC]
* - Lapse rate "t_lapse_rate" [degC/m]
*
* Provides:
* - Air temperature "t" [degC]
* - Lapse rate "t_lapse_rate" [degC/m]
*
*/
class Dist_tlapse : public module_base
{
public:
    Dist_tlapse(config_file cfg);
    ~Dist_tlapse();
    virtual void run(mesh_elem& face, boost::shared_ptr<global> global_param);
    virtual void init(mesh domain, boost::shared_ptr<global> global_param);
    struct data : public face_info
    {
        interpolation interp;
    };
};

/**
@}
*/
