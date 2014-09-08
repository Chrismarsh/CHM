#pragma once


#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include "TPSpline.hpp"

#include <cstdlib>
#include <string>

#include <cmath>

#include <math.h>

/**
* \class rh_llra_var
* \brief Linear lapse rate adjust for relative humidity
*
* Monthly-variable linear lapse rate adjustment for relative humidity based upon Liston, et al. (2006).
*
* Depends:
* -  Air temperature "t" [deg C]
*
* Provides:
* - Relative humidity "rh" [%]
* - Relative humidity "rh_llra_var" [%]
*
* Reference:
* > Liston, Glen E., and Kelly Elder. 2006. A meteorological distribution system for high-resolution terrestrial modeling (MicroMet). Journal of hydrometeorology 7: 217-234
*
* Month | Lapse (m^-1)
* ------|-------
* Jan | 0.00041
* Feb | 0.00042
* Mar | 0.00040
* Apr | 0.00039
* May | 0.00038
* Jun | 0.00036
* Jul | 0.00033
* Aug | 0.00033
* Sep | 0.00036
* Oct | 0.00037
* Nov | 0.00040
* Dec | 0.00040

*/
class rh_llra_var : public module_base
{
public:
    rh_llra_var(std::string ID);
    ~rh_llra_var();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);


};


