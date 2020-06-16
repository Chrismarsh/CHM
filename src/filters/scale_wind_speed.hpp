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

#include "filter_base.hpp"
#include <physics/Atmosphere.h>


/**
 * \ingroup filters wind
 * @{
 * \class scale_wind_speed
 * Scales a wind speed from a given height to the standard reference height in CHM, 50 m
 * Uses a snow roughness z0, but assumes no snow depth (can be corrected for later).
 *
 *
 * **Requires:**
 * - Windspeed [\f$m \cdot s^{-1} \f$ ]
 *
 * **Modifies:**
 * - None
 *
 * **Provides:**
 * - Reference windspeed at reference height  - "U_R"  [\f$m \cdot s^{-1} \f$ ]
 *
 * **Configuration keys:**
 * \rst
 * .. code:: json
 *
 *    {
 *        "variable":"u",
 *        "Z_F":10
 *    }
 *
 * .. confval:: variable
 *
 *    :type: string
 *
 *    Name of the wind variable in the met file
 *
 * .. confval:: Z_F
 *
 *    :type: double
 *    :units: [:math:`m`]
 *
 *    Wind observation height
 *
 * \endrst
 * @}
 */
class scale_wind_speed : public filter_base
{
REGISTER_FILTER_HPP(scale_wind_speed);
private:

    double Z_F; // Measurement height
    double Z_R; // Reference height
    std::string var; //variable name
public:
    scale_wind_speed(config_file cfg);
    ~scale_wind_speed();
    void init();
    void process(std::shared_ptr<station>& station);
};
