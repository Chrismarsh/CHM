#pragma once

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

/**
* \addtogroup modules
* @{
* \class rh_no_lapse
* \brief RH interpolation
*
*
* Depends from met:
* - Relative Humidity 'rh' [%]
*
* Provides:
* - Relative Humidity 'rh' [%]
*/
class rh_no_lapse : public module_base
{
public:
    rh_no_lapse(config_file cfg);

    ~rh_no_lapse();

    virtual void run(mesh_elem &face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};



/**
@}
*/
