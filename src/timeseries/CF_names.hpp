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

#include <string>
#include <boost/bimap.hpp>
#include <boost/assign.hpp>

// This maps netCDF CF names to/from internal CHM names
// The bimap does not support duplicate keys so only supports the current mapping
namespace CF_name_mapping
{
  // two way bimap association
    typedef boost::bimap<std::string, std::string> _t_cf_bimap;
    inline const _t_cf_bimap standard_names = boost::assign::list_of<_t_cf_bimap::relation>
          ("t", "air_temperature")
          ("rh", "relative_humidity")
          ("t_lapse_rate","air_temperature_lapse_rate")
          ("vw_dir", "wind_from_direction")
          ("U_R", "wind_speed")
          ("press", "surface_air_pressure")
          // ("Qli", "downwelling_longwave_flux_in_air")
          // ("Qli", "downwelling_longwave_flux")
          // ("Qli", "surface_downwelling_longwave_flux_in_air")
          ("Qli", "surface_downwelling_longwave_flux")
          ("Qsi", "surface_downwelling_shortwave_flux")
          // ("Qsi", "surface_downwelling_shortwave_flux_in_air")
          // ("Qsi", "downwelling_shortwave_flux")
          // ("Qsi", "downwelling_shortwave_flux_in_air")
          ("p", "precipitation_amount");


}