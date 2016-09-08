#pragma once

#include "filter_base.h"
#include <../modules/constants/Atmosphere.h>


/**
 * \addtogroup filters
 * @{
 * \class scale_wind_speed
 * \brief Scales station/model grid cell wind speed from measured/modeled height to standard reference height for CHM
 *
 * Example call in config file
 * "filter":
       {
         "scale_wind_speed":
         {
           "variable":"u",
           "Z_R":50  // [m]
         }
       }
 *
 * Depends:
 * - U_F [m/s]
 * - Z_F [m] - Height of wind speed measurement/model layer
 *
 */
class scale_wind_speed : public filter_base
{

public:
    scale_wind_speed();
    ~scale_wind_speed();
    void process(boost::shared_ptr<station> station);
};
