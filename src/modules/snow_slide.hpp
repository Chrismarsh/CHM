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
class snow_slide : public module_base
{
public:
    snow_slide(config_file cfg);

    ~snow_slide();

    virtual void run(mesh domain);

    virtual void init(mesh domain);

    void checkpoint(mesh domain,  netcdf& chkpt);
    void load_checkpoint(mesh domain,  netcdf& chkpt);

    struct data : public face_info
    {
        double maxDepth; // m
        double snowdepthavg_copy; // m
        double swe_copy; // m (Note: swe units outside of snowslide are still mm)
        double delta_avalanche_snowdepth; // m^3
        double delta_avalanche_mass; // m^3
    };

};



