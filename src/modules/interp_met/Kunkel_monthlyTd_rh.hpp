#pragma once


#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include <cstdlib>
#include <string>
#include <cmath>
#include <math.h>

#include <meteoio/MeteoIO.h>

/**
* \addtogroup modules
* @{
* \class Kunkel_monthlyTd_rh
* \brief Linear lapse rate adjust for relative humidity
*
* Monthly-variable linear lapse rate adjustment for relative humidity based upon Kunkel. RH is lapsed via dew point temperatures
*
* Depends:
* -  Air temperature "t" [deg C]
*
 * Depends from Met:
 * - Relative Humidity "rh" [%]
 *
* Provides:
* - Relative humidity "rh" [%]
*
* Reference:
* Kunkel, K. E. (1989). Simple procedures for extrapolation of humidity variables in the mountainous western United States. Journal of Climate, 2(7), 656–669. Retrieved from http://ams.allenpress.com/perlserv/?request=get-abstract&amp;doi=10.1175/1520-0442(1989)002<0656:SPFEOH>2.0.CO;2
*/
class Kunkel_monthlyTd_rh : public module_base
{
public:
    Kunkel_monthlyTd_rh(config_file cfg);
    ~Kunkel_monthlyTd_rh();
    virtual void run(mesh_elem& elem, boost::shared_ptr<global> global_param);
    virtual void init(mesh domain, boost::shared_ptr<global> global_param);
    struct data : public face_info
    {
        interpolation interp;
    };
};


/**
@}
*/