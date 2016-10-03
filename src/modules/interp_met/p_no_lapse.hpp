#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include <cstdlib>
#include <string>

#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>
#include <gsl/gsl_fit.h>
/**
* \addtogroup modules
* @{
* \class Precip
* \brief Calculates precipitation
*
* Spatially distributes liquid water precipitation using a constant lapse rate
*
* Depends:
* - Precip from met file "p" [mm]
*
* Provides:
* - Precip "p" [mm]
*/
class p_no_lapse : public module_base
{
public:
    p_no_lapse(config_file cfg);
    ~p_no_lapse();
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