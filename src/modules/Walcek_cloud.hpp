#pragma once

#include "module_base.hpp"
#include <math.h>
#include <algorithm>


#include <meteoio/MeteoIO.h>
/**
* \addtogroup modules
* @{
* \class Walcek_cloud
* \brief Calculates cloud fraction
*
* Calculates a cloud fraction
* Walcek, C. J. (1994). Cloud cover and its relationship to relative humidity during a springtime midlatitude cyclone. Monthly Weather Review, 122(6), 1021–1035.
*
* Depends:
* - Air temperature (t)
* - Relative humidity (rh)
*
* Provides:
* - Atmospheric transmittance "cloud_frac" [-]
*/
class Walcek_cloud : public module_base
{
public:
    Walcek_cloud();
    ~Walcek_cloud();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);

};

/**
@}
*/