//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once

#include <boost/shared_ptr.hpp>

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"


#include <string>

/**
 * \ingroup modules snow precipitation
 * @{
 * \class threshold_p_phase
 *
 * Calculates phase from air temperature. Defaults to 2 \f$ {}^\circ C \f$, looks for the config parameter "threshold_temperature".
 * Binary snow/rain at 2 \f$ {}^\circ C \f$ threshold.
 *
 * **Depends:**
 * - Air temperature "t" [ \f$ {}^\circ C \f$ ]
 * - Precipitation "p" [ \f$ mm \cdot dt^{-1} \f$ ]
 *
 * **Provides:**
 * - Fractional precipitation that is rain "frac_precip_rain" [0-1]
 * - Fractional precipitation that is snow "frac_precip_snow" [0-1]
 * - Precipitation, rain "p_rain" [\f$ mm \cdot dt^{-1} \f$]
 * - Precipitation, snow "p_snow" [\f$ mm \cdot dt^{-1} \f$]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "threshold_temperature": 2.0
 *    }
 *
 * .. confval:: threshold_temperature
 *
 *    :default: 2.0  :math:`{}^\circ C`
 *
 *    Threshold air temperature for rain/snow delineation
 *
 * \endrst
 *
 * @}
 */
class threshold_p_phase : public module_base
{
REGISTER_MODULE_HPP(threshold_p_phase);
public:
    threshold_p_phase(config_file cfg);

    ~threshold_p_phase();

    virtual void run(mesh_elem &face);

private:

    // the air temperature threshold above which the precip phase is liquid
    double t_thresh;

};
