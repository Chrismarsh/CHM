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
#include <math.h>


/**
 * \ingroup filters precip
 * @{
 * \class macdonald_undercatch
 * Undercatch correction for a Alter shielded Geonor and tipping bucket via Macdonald, et al. 2007
 *
 * **Requires:**
 * - Precipitation [\f$ mm\f$ ]
 * - Windspeed at gauge height  [\f$m \cdot s^{-1} \f$ ]
 *
 * **Modifies**:
 * - Input precipitation
 *
 * **Configuration keys:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "precip_var": "ppt",
 *       "wind_var": "u"
 *    }
 *
 * .. confval:: precip_var
 *
 *    :type: string
 *
 *    Name of the precipitation variable in the met file.
 *
 * .. confval:: wind_var
 *
 *    :type: string
 *
 *    Name of the wind speed variable in the met file.
 *
 * \endrst
 *
 * **References:**
 * - Macdonald, J., & Pomeroy, J. (2007). Gauge Undercatch of Two Common Snowfall Gauges in a Prairie Environment. Proceedings of the 64th Eastern Snow Conference, St. John‘s, Canada., 119–126.
 * @}
 * */
class macdonald_undercatch : public filter_base
{
REGISTER_FILTER_HPP(macdonald_undercatch);
private:
    std::string precip_var;
    std::string wind_var;
public:
    macdonald_undercatch(config_file cfg);
    ~macdonald_undercatch();
    void init();
    void process(std::shared_ptr<station>& station);
};
