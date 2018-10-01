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


#include "fast_shadow.hpp"
REGISTER_MODULE_CPP(fast_shadow);

fast_shadow::fast_shadow(config_file cfg)
        : module_base(parallel::data)
{
    depends("solar_az");
    depends("solar_el");
    provides("shadow");

    //number of steps along the search vector to check for a higher point
    steps = cfg.get("steps",10);
    //max distance to search
    max_distance = cfg.get("max_distance",1000.0);

    //size of the step to take
    size_of_step = max_distance / steps;
}

fast_shadow::~fast_shadow()
{


}

void fast_shadow::run(mesh_elem& face)
{

    double solar_el = face->face_data("solar_el") *M_PI / 180.;

    face->set_face_data("shadow", 0);
    //bail early
    if (solar_el < 0)
        return;

    double solar_az = face->face_data("solar_az") ;

    Point_3 me = face->center();

    double phi = 0.;
    // search along each azimuth in j step increments to find horizon angle
    for (int j = 1; j <= steps; ++j)
    {
        double distance = j * size_of_step;

        auto f = face->find_closest_face(solar_az, distance);

        double z_diff = f->center().z() - me.z() ;
        if (z_diff > 0)
        {
            double dist = math::gis::distance(f->center(), me);
            phi = std::max(atan(z_diff / dist), phi);
        }
        //try to bail early if possible
        if (phi > solar_el )
        {
            face->set_face_data("shadow", 1);
        }
    }
    //try to bail early if possible
    if (phi > solar_el )
    {
        face->set_face_data("shadow", 1);
    }



}
