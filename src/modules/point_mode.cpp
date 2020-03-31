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

#include "point_mode.hpp"
REGISTER_MODULE_CPP(point_mode);

point_mode::point_mode(config_file cfg)
        : module_base("point_mode", parallel::data, cfg)
{

     t              = cfg.get("provide.t",true);
     rh             = cfg.get("provide.rh",true);
     U_R            = cfg.get("provide.U_R",true);
     U_2m_above_srf = cfg.get("provide.U_2m_above_srf",true);
     p              = cfg.get("provide.p",true);
     ilwr           = cfg.get("provide.ilwr",true);
     iswr           = cfg.get("provide.iswr",true);
     vw_dir         = cfg.get("provide.vw_dir",true);
     iswr_diffuse   = cfg.get("provide.iswr_diffuse",false);
     iswr_direct    = cfg.get("provide.iswr_direct",false);
     T_g    = cfg.get("provide.T_g",false);

    if(t)
    {
        depends_from_met("t");
        provides("t");
    }

    if(rh)
    {
        depends_from_met("rh");
        provides("rh");
    }

    // If U_2m_above_srf is provided, use it. Otherwise use U_R.
    if(U_2m_above_srf) {
        depends_from_met("u");
        provides("U_2m_above_srf");
    }else if(U_R) {
        depends_from_met("U_R");
        provides("U_R");
    }

    if(vw_dir)
    {
        depends_from_met("vw_dir");
        provides("vw_dir");
    }

    if(p)
    {
        depends_from_met("p");
        provides("p");
    }

    if(ilwr)
    {
        depends_from_met("Qli");
        provides("ilwr");
    }

    if(iswr)
    {
        depends_from_met("Qsi");
        provides("iswr");
    }

    if(iswr_diffuse)
    {
        depends_from_met("iswr_diffuse");
        provides("iswr_diffuse");
    }

    if(iswr_direct)
    {
        depends_from_met("iswr_direct");
        provides("iswr_direct");
    }

    if(T_g)
    {
        depends_from_met("T_g");
        provides("T_g");
    }

}

point_mode::~point_mode()
{

}

void point_mode::run(mesh_elem &face)
{
    // at this point, if the user has provided more than 1 station, they've been stripped out.
    // we can safetly take the 1st (and only) station.

    if(t)
    {
        double st =(*face->nearest_station())["t"_s];
        (*face)["t"_s]=st;
    }

    if(rh)
    {
        double srh =(*face->nearest_station())["rh"_s];
        (*face)["rh"_s]= srh;
    }

    if(U_2m_above_srf) {
        double su =(*face->nearest_station())["u"_s];
        su = std::max(su,0.1);
        (*face)["U_2m_above_srf"_s]=su;
    } else if (U_R)
    {
        double su =(*face->nearest_station())["U_R"_s];

        //make sure we don't have zero wind speeds
        su = std::max(su,0.1);
        (*face)["U_R"_s]=su;
    }

    if(vw_dir)
    {
        double sdir =(*face->nearest_station())["vw_dir"_s];
        (*face)["vw_dir"_s]=sdir;
    }

    if(p)
    {
        double sp =(*face->nearest_station())["p"_s];
        (*face)["p"_s]= sp;
    }
    if(ilwr)
    {
        double silwr =(*face->nearest_station())["Qli"_s];
        (*face)["ilwr"_s]= silwr;
    }
    if(iswr)
    {
        double iswr =(*face->nearest_station())["Qsi"_s];
        (*face)["iswr"_s]= iswr;

    }
    if(iswr_diffuse)
    {
        double iswr_diffuse =(*face->nearest_station())["iswr_diffuse"_s];
        (*face)["iswr_diffuse"_s]= iswr_diffuse;

    }
    if(iswr_direct)
    {
        double iswr_direct =(*face->nearest_station())["iswr_direct"_s];
        (*face)["iswr_direct"_s]= iswr_direct;

    }

    if(T_g)
    {
        double T_g =(*face->nearest_station())["T_g"_s];
        (*face)["T_g"_s]= T_g;
    }

    face->parameter("svf") = cfg.get("override.svf", face->parameter("svf"_s));

}
