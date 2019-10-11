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

#include "scale_wind_speed.hpp"
REGISTER_FILTER_CPP(scale_wind_speed);

scale_wind_speed::scale_wind_speed(config_file cfg)
  : filter_base("scale_wind_speed", cfg)
{

}

scale_wind_speed::~scale_wind_speed()
{



}

void scale_wind_speed::init(std::shared_ptr<station>& station)
{
    //look at the config data to determine what we are modifying
    var = cfg.get<std::string>("variable");
    Z_F          = cfg.get<double>("Z_F"); // Get measurement height [m]
    Z_R          = Atmosphere::Z_U_R; // Reference wind speed height [m]

    // Initialize new wind speed at ref height variable
    station->add_variable("U_R");
}
void scale_wind_speed::process(std::shared_ptr<station>& station)
{
    double U_F = station->now().get(var); // Here wind u [m/s] at Z_U
    double U_R = -9999;
    if(!is_nan(U_F))
    {
        U_R = Atmosphere::log_scale_wind(U_F, Z_F, Z_R, 0); // Assume 0 snow depth
    }

    station->now().set("U_R",U_R);

}
