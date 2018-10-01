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

#include "sub_grid.hpp"
REGISTER_MODULE_CPP(sub_grid);

sub_grid::sub_grid(config_file cfg)
        :module_base(parallel::data)

{
    depends("snowdepthavg");

    provides("snowcoverfraction");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}

sub_grid::~sub_grid()
{


}

void sub_grid::run(mesh_elem &face) {

    double snowdepthavg = face->face_data("snowdepthavg");
    double snowcoverfraction;

    if(snowdepthavg>0) {
        snowcoverfraction = 1;
    } else {
        snowcoverfraction = 0;
    }

    face->set_face_data("snowcoverfraction",snowcoverfraction);

}
