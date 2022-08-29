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

#include "mpi_example.hpp"
REGISTER_MODULE_CPP(mpi);

mpi::mpi(config_file cfg)
        : module_base("mpi", parallel::domain, cfg)
{
    provides("mpi_rank");
    provides("mpi_ghost_test");
}

mpi::~mpi()
{

}

void mpi::run(mesh& domain)
{

    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i); // Get face
        (*face)["mpi_rank"] = face->owner;
        for (int j = 0; j < 3; ++j)
        {
            auto n = face->neighbor(j);
            if (n != nullptr && n->is_ghost)
            {
                (*n)["mpi_ghost_test"] = face->owner;
            }
        }
    }

    // This will set the ghosts on our faces to equal the mpi rank they came from
    // their face -> our ghost
    domain->ghost_neighbors_communicate_variable("mpi_rank");

    // this will set the faces on the remote ranks to be our rank
    // our ghost -> their face
    domain->ghost_to_neighbors_communicate_variable("mpi_ghost_test");
}
