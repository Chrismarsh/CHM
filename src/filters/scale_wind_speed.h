/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "filter_base.h"
#include <constants/Atmosphere.h>


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
private:
    double Z_F;
    double Z_R;
    std::string var;
public:
    scale_wind_speed();
    ~scale_wind_speed();
    void init(boost::shared_ptr<station>& station);
    void process(boost::shared_ptr<station>& station);
};
