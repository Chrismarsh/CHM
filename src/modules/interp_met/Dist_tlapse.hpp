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
* \class Liston_monthly_llra_ta
* \brief Constant monthly linear lapse rate adjustment
*
* Constant monthly linear lapse rate adjustment for air temperature.
*
* Month | Lapse rate (degC/m)
* ------|----------------
* Jan | 0.0044
* Feb | 0.0059
* Mar | 0.0071
* Apr | 0.0078
* May | 0.0081
* Jun | 0.0082
* Jul | 0.0081
* Aug | 0.0081
* Sep | 0.0077
* Oct | 0.0068
* Nov | 0.0055
* Dec | 0.0047
*
* Depends:
* - None
*
* Provides:
* - Air temperature "t" [degC]
* - Air temperatue "Liston_monthly_llra_ta" [degC]
*
* Reference:
* > Liston, Glen E., and Kelly Elder. 2006. A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). Journal of hydrometeorology 7: 217-234

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
