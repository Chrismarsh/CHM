#pragma once

#include "module_base.hpp"
#include <math.h>
#include <algorithm>


#include <meteoio/MeteoIO.h>
/**
* \addtogroup modules
* @{
* \class Walcek_atm_trans
* \brief Calculates shortwave atmospheric transmittance
*
* Calculates incoming direct-beam shortwave solar radiation transmittance following
* Walcek, C. J. (1994). Cloud cover and its relationship to relative humidity during a springtime midlatitude cyclone. Monthly Weather Review, 122(6), 1021–1035.
*
* Depends:
* - Air temperature (t)
* - Relative humidity (rh)
*
* Provides:
* - Atmospheric transmittance "atm_trans" [-]
*/
class Walcek_atm_trans : public module_base
{
public:
    Walcek_atm_trans();
    ~Walcek_atm_trans();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);

};

/**
@}
*/