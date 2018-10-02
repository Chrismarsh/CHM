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

#include "iswr.hpp"
REGISTER_MODULE_CPP(iswr);

iswr::iswr(config_file cfg)
        :module_base("iswr", parallel::data, cfg)
{

    depends("iswr_diffuse_no_slope");
    depends("iswr_direct_no_slope");

    depends("solar_el");
    depends("solar_az");


    provides("iswr");
    provides("iswr_direct");
    provides("iswr_diffuse");

    provides("solar_angle");

    optional("shadow");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;

    assume_no_slope = cfg.get("no_slope",false);
    already_cosine_corrected = cfg.get("already_cosine_corrected",false);


}
void iswr::run(mesh_elem& face)
{

    if(global_param->is_point_mode() && !already_cosine_corrected)
    {
        LOG_ERROR << "Most observations implicitly have a cosine-correction built in by virtu of the flat-plane observation. "
                     "When using point-mode, you probably want to set -c config.iswr.already_cosine_corrected:true";
    }

    double A = face->face_data("solar_az") * mio::Cst::to_rad;
    double E = face->face_data("solar_el") * mio::Cst::to_rad;

    //radiation data
    //solar vector
    //xyz cartesian
    arma::vec S;
    S << cos(E) * sin(A) << arma::endr
      << cos(E) * cos(A) << arma::endr
      << sin(E) << arma::endr;

    arma::vec N;
    if (assume_no_slope)
    {

        N << 0 << arma::endr
        << 0 << arma::endr
        << 1 << arma::endr;
    }
    else
    {
        Vector_3 n = face->normal();


        N << n[0] << arma::endr
        << n[1] << arma::endr
        << n[2] << arma::endr;
    }


    double angle = acos(arma::dot(S,N));
    angle = cos(angle);

    if(angle < 0.0 || E < 0.0523598776) //3deg -> rad
        angle = 0.0;

    face->set_face_data("solar_angle",angle);

    //if we have remote shadowing
    if(has_optional("shadow"))
    {
        if(face->face_data("shadow") == 1)
        {
            angle = 0;
        }
    }

    double direct_beam = face->face_data("iswr_direct_no_slope");

    //If we're using obs at a point, this should be set to true
    if(already_cosine_corrected)
    {
        direct_beam = direct_beam / sin(E);
    }

    double swr =  angle * direct_beam;
    double diff = face->face_data("iswr_diffuse_no_slope");

    swr = std::max(0.0, swr);
    diff = std::max(0.0,diff);

    face->set_face_data("iswr_direct",swr );
    face->set_face_data("iswr_diffuse",diff );
    face->set_face_data("iswr", swr + diff );

}

iswr::~iswr()
{


}
