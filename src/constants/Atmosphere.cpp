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

#include <constants/Atmosphere.h>

namespace Atmosphere
{

    // Logrithmic, assuming no snow cover and no canopy bewteen Z_in and Z_out
    // Because filter is before runtime, we do not know what the snowdepth will be, thus this
    // introduces some error latter when wind is scaled down taking into account the snowdepth.
    double log_scale_wind(double u, double Z_in, double Z_out, double snowdepthavg, double z0)
    {
        double Z0_SNOW  = z0; //Snow::Z0_SNOW; // Snow roughness (m)

        u = u * log((Z_out - (snowdepthavg + Z0_SNOW)) / Z0_SNOW) / log((Z_in - (snowdepthavg + Z0_SNOW)) / Z0_SNOW);
        return u;
    }

    // Exponential following Inoue (1963)
    double exp_scale_wind(double u, double Z_in, double Z_out, const double alpha)
    {

        u = u * exp(alpha*(Z_out/Z_in-1));
        return u;
    }

    // Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
    // See Fig 1 and 2 of Kienzle (2010, Hydrological Processes)  
    // Slope in radian
    double corr_precip_slope(double p, double slope)
    {
        p = p * cos(slope);
        return p;
    }


}
