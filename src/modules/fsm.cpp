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

#include "fsm.hpp"
REGISTER_MODULE_CPP(FSM);

FSM::FSM(config_file cfg)
    : module_base("FSM", parallel::data, cfg)
{
    depends("solar_el");
    depends("ilwr");
    depends("rh");
    depends("t");
    depends("p_snow");
    depends("p_rain");
    depends("U_2m_above_srf");

    depends("iswr_direct");
    depends("iswr_diffuse");

    optional("iswr_subcanopy");
    optional("rh_subcanopy");
    optional("ta_subcanopy");
    optional("ilwr_subcanopy");

    provides("swe");
    provides("snowdepthavg");
    provides("snowdepthavg_vert");

    conflicts("snow_slide"); // for now don't run w/ snowslide
}

void FSM::init(mesh& domain)
{
    __layers_MOD_ncnpy = 2;
    __layers_MOD_nsmax = 3;
    __layers_MOD_nsoil = 4;

    allocate();

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<data>(ID);

    }
}
void FSM::run(mesh_elem& face)
{
    if(is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }

    // met data

    float zT = 2; // m
    float zU = 2;
    float Ps = mio::Atmosphere::stdAirPressure(face->get_z()); // Pa

    float dt = (float)global_param->dt();


    float rh = -9999;
    if(has_optional("rh_subcanopy")) {
        rh = (*face)["rh_subcanopy"_s];
    } else {
        rh = (*face)["rh"_s];
    }

    //TODO: needs to have subcanopy added?
    float Rf = (*face)["p_rain"_s] / dt; // rainfall rate
    float Sf = (*face)["p_snow"_s] / dt; // snowfall rate

    auto d = face->get_module_data<data>(ID);

    float elev = (float)(*face)["solar_el"_s];
    float ilwr = -9999;
    if(has_optional("ilwr_subcanopy")) {
        ilwr = (float)(*face)["ilwr_subcanopy"_s];
    } else {
        ilwr = (float)(*face)["ilwr"_s];
    }


    // TODO: needs subcanopy
    float Sdiff = (float)(*face)["iswr_diffuse"_s];
    float Sdir = (float)(*face)["iswr_direct"_s];

    //TODO: needs avalanching mass

    float t = -9999;
    if(has_optional("ta_subcanopy")) {
        t = (float)(*face)["ta_subcanopy"_s];
    } else {
        t = (float)(*face)["t"_s];
    }
    t += 271.15;

    float tc = t - 273.15;
    float Qs = __constants_MOD_eps * (__constants_MOD_e0 / Ps) *
               exp((float)17.5043 * tc / ((float)241.3 + tc));
    float Qa = (rh/ (float)100.0) * Qs; // specific humidity

    float U = (float)(*face)["U_2m_above_srf"_s];


    fsm2_timestep(
        // Driving variables
        &dt, &elev, &zT, &zU, &ilwr,
        &Ps, &Qa, &Rf, &Sdiff, &Sdir, &Sf, &t, &U,

        // Vegetation characteristics
        &d->veg.alb0, &d->veg.hveg, &d->veg.VAI,

        // State variables
        &d->state.albs, &d->state.Tsrf, d->state.Dsnw, &d->state.Nsnow, d->state.Qcan,
        d->state.Rgrn, d->state.Sice, d->state.Sliq, d->state.Sveg, d->state.Tcan, d->state.Tsnow,
        d->state.Tsoil, d->state.Tveg, d->state.Vsmc,

        // Diagnostics
        &d->diag.H, &d->diag.LE, &d->diag.LWout, &d->diag.LWsub, &d->diag.Melt,
        &d->diag.Roff, &d->diag.snd, &d->diag.Svg, &d->diag.SWE, &d->diag.SWout, &d->diag.SWsub,
        &d->diag.Usub);

    (*face)["swe"_s] = d->diag.SWE;
    (*face)["snowdepthavg"_s] = d->diag.snd;
    (*face)["snowdepthavg_vert"_s] = d->diag.snd/std::max(0.001,cos(face->slope()));

}


FSM::~FSM()
{


}
