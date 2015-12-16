#pragma once

#include "module_base.hpp"
#include <gsl/gsl_fit.h>
#include <vector>
#include <gsl/gsl_combination.h>

/**
* \addtogroup modules
* @{
* \class Precip
* \brief Calculates precipitation
*
* Spatially distributes liquid water precipitation using the Thornton, et al. 1997 method.
* Unlike Thornton_p, this calculates scaling rates based on observed data on a per-timestep basis
*
*
* Depends:
* - Precip from met file "p" [mm]
*
* Provides:
* - Precip "p" [mm]
*
* Reference:
* - Thornton, P. E., Running, S. W., & White, M. A. (1997). Generating surfaces of daily meteorological variables over large regions of complex terrain. Journal of Hydrology, 190(3-4), 214–251. http://doi.org/10.1016/S0022-1694(96)03128-9
* */
class Thornton_var_p : public module_base
{
public:
    Thornton_var_p(config_file cfg);
    ~Thornton_var_p();
    void run(mesh_elem& elem, boost::shared_ptr<global> global_param);
    virtual void init(mesh domain, boost::shared_ptr<global> global_param);
    struct data : public face_info
    {
        interpolation interp;
    };
};

