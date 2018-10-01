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

#include "threshold_p_phase.hpp"
REGISTER_MODULE_CPP(threshold_p_phase);

threshold_p_phase::threshold_p_phase(config_file cfg)
    : module_base(parallel::data)
{
    depends("t");
    depends("p");

    t_thresh = cfg.get("threshold_temperature", 2.0) ;

    provides("frac_precip_rain");
    provides("frac_precip_snow");

    provides("p_snow");
    provides("p_rain");

}

threshold_p_phase::~threshold_p_phase()
{

}

void threshold_p_phase::run(mesh_elem &face)
{
    double p = face->face_data("p");

    if( face->face_data("t") >= t_thresh)
    {
        face->set_face_data("frac_precip_rain", 1);
        face->set_face_data("frac_precip_snow", 0);

        face->set_face_data("p_rain", p);
        face->set_face_data("p_snow", 0 );
    } else
    {
        face->set_face_data("frac_precip_rain", 0);
        face->set_face_data("frac_precip_snow", 1);

        face->set_face_data("p_rain", 0);
        face->set_face_data("p_snow", p );
    }



}
