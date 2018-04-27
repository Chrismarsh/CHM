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

/**
* \addtogroup modules
* @{
* \class kunkel_rh
* \brief RH interpolation
*
* Interpolates RH with elevation. Uses eqn 15 in Kunkel to directly interplate RH.
*
* Depends from met:
* - Relative Humidity 'rh' [%]
*
* Provides:
* - Relative Humidity 'rh' [%]
 *
 * References:
 * - Kunkel, K. E. (1989). Simple procedures for extrapolation of humidity variables in the mountainous western United States. Journal of Climate, 2(7), 656–669. Retrieved from http://ams.allenpress.com/perlserv/?request=get-abstract&amp;doi=10.1175/1520-0442(1989)002<0656:SPFEOH>2.0.CO;2
*/
class kunkel_rh : public module_base
{
public:
    kunkel_rh(config_file cfg);

    ~kunkel_rh();

    virtual void run(mesh_elem &face);
    virtual void init(mesh domain);
    struct data : public face_info
    {
        interpolation interp;
    };
};



/**
@}
*/