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

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"

#include <cstdlib>
#include <string>
#include <cmath>
#include <armadillo>
#define _USE_MATH_DEFINES
#include <math.h>

#include <boost/math/tools/roots.hpp>


/**
 * \ingroup modules precip snow
 * @{
 * \class Harder_precip_phase
 * Calculates precipitation phase via falling hydrometeor energy balance following Harder, et al (2013)
 *
 * **Depends:**
 * - Air temperature "t" [ \f$ {}^\circ C\f$ ]
 * - Relative Humidity "rh"  [%]
 * - Precip "p" [ \f$mm \cdot dt^{-1}\f$ ]
 *
 * **Provides:**
 * - Snow precip p_snow [\f$mm \cdot dt^{-1}\f$]
 * - Liquid precip p_rain [\f$mm \cdot dt^{-1}\f$]
 * - Fraction of rain "frac_precip_rain" [-]
 * - Fraction of snow "frac_precip_snow" [-]
 * - Cumulated snow precip "acc_snow" [mm]
 * - Cumulated liquid precip "acc_rain" [mm]
 *
 *
 * **References:**
 * - Harder, P., Pomeroy, J. (2013). Estimating precipitation phase using a psychrometric energy balance method
 * Hydrological Processes  27(13), 1901-1914. https://dx.doi.org/10.1002/hyp.9799
 *
 * @}
*/
class Harder_precip_phase : public module_base
{
REGISTER_MODULE_HPP(Harder_precip_phase);
public:
    Harder_precip_phase(config_file cfg);
    ~Harder_precip_phase();
    virtual void run(mesh_elem& face);
    void init(mesh& domain);
    double b;
    double c;

    class data : public face_info
    {
    public:
        double hours_since_snowfall;
        double acc_rain;
        double acc_snow;

    };

    void checkpoint(mesh& domain,  netcdf& chkpt);
    void load_checkpoint(mesh& domain, netcdf& chkpt);


};

