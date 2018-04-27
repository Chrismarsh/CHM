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

#include "slope_iswr.hpp"

slope_iswr::slope_iswr(config_file cfg)
        :module_base(parallel::data)
{
    depends("iswr");
    depends("iswr_diffuse");
    depends("iswr_direct");

    depends("solar_el");
    depends("solar_az");


    provides("iswr");
    provides("iswr_direct");
    provides("solar_angle");

    optional("shadow");

    LOG_DEBUG << "Successfully instantiated module " << this->ID;

    assume_no_slope = cfg.get("no_slope",false);
}
void slope_iswr::run(mesh_elem& face)
{
    

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

    double swr =  angle * face->face_data("iswr_direct");
    double diff = face->face_data("iswr_diffuse");

    face->set_face_data("iswr_direct",swr );
    face->set_face_data("iswr", swr + diff );

}

slope_iswr::~slope_iswr()
{
    
    
}



