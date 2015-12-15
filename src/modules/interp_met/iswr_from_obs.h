#pragma once
#include "module_base.hpp"
#include <math.h>
#include <algorithm>
#include <meteoio/MeteoIO.h>


/**
* \addtogroup modules
* @{
* \class iswr_from_obs
* \brief Spatially interpolates observed shortwave measurements
*
* Spatially interpolates observed shortwave measurements and splits it into direct and diffuse beams. This split
 * is via Iqbal
*
* Depends:
* -  Shortwave radiation met file "Qsi" [W/m^2]
*
* Provides:
* - Shortwave all beam "iswr" [W/m^2]
* - Shortwave direct "iswr_direct" [W/m^2]
* - Shortwave diffuse "iswr_diffuse" [W/m^2]
* - Atmospheric transmittance, [0,1] "atm_trans" [-]
*/
class iswr_from_obs : public module_base
{
public:
    iswr_from_obs();
    ~iswr_from_obs();
    void run(mesh_elem &elem, boost::shared_ptr<global> global_param);
    virtual void init(mesh domain, boost::shared_ptr<global> global_param);
    struct data : public face_info
    {
        interpolation interp;
    };
};


