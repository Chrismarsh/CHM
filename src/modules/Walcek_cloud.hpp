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
* Calculates a cloud fraction. Extrapolates RH700mb via Kunkel
* Walcek, C. J. (1994). Cloud cover and its relationship to relative humidity during a springtime midlatitude cyclone. Monthly Weather Review, 122(6), 1021–1035.
* Kunkel, K. E. (1989). Simple procedures for extrapolation of humidity variables in the mountainous western United States. Journal of Climate, 2(7), 656–669. Retrieved from http://ams.allenpress.com/perlserv/?request=get-abstract&amp;doi=10.1175/1520-0442(1989)002<0656:SPFEOH>2.0.CO;2
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