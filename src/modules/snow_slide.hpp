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
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_sort.h>
#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"

#include <string>

/**
 * \ingroup modules snow
 * @{
 * \class snow_slide
 *
 * SnowSlide is a simple topographically-driven model that simulates the effects of gravitational snow transport.
 * It uses a snow holding depth that decreases exponentially with increasing slope angle, limiting snow accumulation
 * in steep terrain. The algorithm moves mass from the highest triangle of the mesh to the lowest one. If the snow
 * depth exceeds the snow holding capacity for a given triangle, excess snow is redistributed to the lower adjacent
 * triangles, proportionally to the elevation difference between the neighboring triangles and the original one.
 * SnowSlide uses the total elevation (snow depth plus surface elevation) to operate. In this study, the default
 * formulation of the snow holding depth proposed by Bernhardt and Schulz (2010) is used which leads to a maximal
 * snow thickness (taken perpendicular to the slope) of 3.08 m, 1.11 m, 0.45 m, and 0.15 m for slopes of 30° 45°, 60°,
 * and 75°, respectively.
 *
 * In the manuscript it is unclear if the holding capacity (max depth) is parameterized for a snow thickness normal to
 * the slope (i.e., how the snowmodels treat snow) or if it is in the vertical direction (i.e., how LiDAR would
 * measure it, a.k.a cosine corrected snow depth). An analysis versus observed LiDAR derived snowdepths over the
 * Kananaskis, Canada domain as well as data observed by Sommer et al. (2015) showed a better evaluation with snowdepths
 * if it is assumed the parameterization is defined for snowdepths in the vertical.
 *
 * \rst
 * .. image:: images/snowslide_eval1.png
 * \endrst
 *
 * **Depends:**
 * - Snow depth "snowdepthavg" [m]
 * - Snow depth "snowdepthavg_vert" [m]
 * - Snow Water Equivalent "swe" [mm]
 *
 * **Provides:**
 * - Change in snow mass due to avalanching "delta_avalanche_mass" [mm]
 * - Change in snow depth due to avalanching "delta_avalanche_snowdepth" [mm]
 * - Maximum snowdepth holding capacity of a triangle "maxDepth" [m]
 *
 * **Configuration:**
 * \rst
 * .. code:: json
 *
 *    {
 *       "avalache_mult": 3178.4,
 *       "avalache_pow": -1.998
 *   }
 *
 * .. confval:: avalache_mult
 *
 *    :default: 3178.4
 *
 * .. confval:: avalache_pow
 *
 *    :default: -1.998
 *
 * \endrst
 *
 * **References:**
 * - Bernhardt, M., Schulz, K. (2010). SnowSlide: A simple routine for calculating gravitational snow transport
 * Geophysical Research Letters  37(11), 1-6. https://dx.doi.org/10.1029/2010gl043086
 * @}
 */
class snow_slide : public module_base
{
REGISTER_MODULE_HPP(snow_slide);
public:
    snow_slide(config_file cfg);

    ~snow_slide();

    virtual void run(mesh& domain);

    virtual void init(mesh& domain);

    void checkpoint(mesh& domain,  netcdf& chkpt);
    void load_checkpoint(mesh& domain,  netcdf& chkpt);

    struct data : public face_info
    {
        double maxDepth; // Vertical snow holding depth  m
        double snowdepthavg_copy; // m
        double snowdepthavg_vert_copy; // m
        double slope; // rad
        double swe_copy; // m (Note: swe units outside of snowslide are still mm)
        double delta_avalanche_snowdepth; // m^3
        double delta_avalanche_mass; // m^3
    };
    bool use_vertical_snow; 
// True: apply the maximal snow holding capacity to snow depth (measured vertically)
// False: apply the maximal snow holding capacity to snow thickness (perpendicular to the surface)

};
