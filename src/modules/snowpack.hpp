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

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"



#include <meteoio/MeteoIO.h>
#include <constants/Atmosphere.h>
#include <snowpack/libsnowpack.h>

#include <string>
class Lehning_snowpack : public module_base
{
REGISTER_MODULE_HPP(Lehning_snowpack);
public:
    Lehning_snowpack(config_file cfg);

    ~Lehning_snowpack();

    virtual void run(mesh_elem &face);

    virtual void init(mesh domain);


    struct data : public face_info
    {
        //main snowpack model
        boost::shared_ptr<Snowpack> sp;

        /*
         * This is the PRIMARY data structure of the SNOWPACK program \n
         * It is used extensively not only during the finite element solution but also to control
         */
        boost::shared_ptr<SnowStation> Xdata;

        boost::shared_ptr<SnowpackConfig> Spackconfig;
        boost::shared_ptr<Meteo> meteo;
        boost::shared_ptr<Stability> stability;
        mio::Config config;
        double cum_precip;

        double sum_subl;
    };

    double sn_dt; // calculation step length
    double const_T_g; // constant ground temp, degC

};
