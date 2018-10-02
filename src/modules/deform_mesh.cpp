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

#include "deform_mesh.hpp"
REGISTER_MODULE_CPP(deform_mesh);

deform_mesh::deform_mesh(config_file cfg)
        : module_base("deform_mesh", parallel::domain, cfg)
{

}

deform_mesh::~deform_mesh()
{

}

void deform_mesh::run(mesh domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_vertex(); i++)
    {
        Point_3 p;

        auto vert = domain->vertex(i);
        double z = vert->point().z();

        if(z > domain->min_z())
        {
           z -= (z-domain->min_z()) * 0.25;
        }


        p = Point_3(vert->point().x(), vert->point().y(), z);
        vert->set_point(p);
    }

    domain->_terrain_deformed = true;
}
