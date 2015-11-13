#pragma once
#include "module_base.hpp"
#include <math.h>
#include <algorithm>
#include <meteoio/MeteoIO.h>
#include "TPSpline.hpp"

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
*/
class iswr_from_obs : public module_base
{
public:
    iswr_from_obs();
    ~iswr_from_obs();
    void run(mesh_elem &elem, boost::shared_ptr<global> global_param);
};


