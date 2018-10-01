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
#include "TPSpline.hpp"

/*
 * When using the model in point mode. Doesn't do any interp, but rather uses a specific input station as input.
 */
class point_mode : public module_base
{
REGISTER_MODULE_HPP(point_mode);
public:
    point_mode(config_file cfg);

    ~point_mode();

    virtual void run(mesh_elem &face);


    bool t ;
    bool rh ;
    bool U_R ;
    bool p ;
    bool ilwr ;
    bool iswr ;
    bool vw_dir ;
    bool iswr_diffuse ;
    bool iswr_direct ;
    bool U_2m_above_srf ;
    bool T_g;


};
