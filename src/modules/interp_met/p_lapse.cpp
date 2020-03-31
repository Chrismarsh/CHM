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

#include "p_lapse.hpp"
REGISTER_MODULE_CPP(p_lapse);

p_lapse::p_lapse(config_file cfg)
        : module_base("p_lapse", parallel::data, cfg)
{

    depends_from_met("p");

    provides("p");
    provides("p_no_slope");


    apply_cosine_correction = cfg.get("apply_cosine_correction",false);


    LOG_DEBUG << "Successfully instantiated module " << this->ID;
}
void p_lapse::init(mesh& domain)
{
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        auto d = face->make_module_data<p_lapse::data>(ID);
        d->interp.init(global_param->interp_algorithm,face->stations().size() );
    }
}
void p_lapse::run(mesh_elem& face)
{

    // Precipitation lapse rate derived from Marmot Creek stations by Logan Fang (used in CRHM for Marmot Creek domain)
    //(100 m)^-1
    double monthly_factors[] = {0.1081,0.1081 ,0.1081, 0.0997, 0.0997, 0.0592, 0.0592, 0.0592, 0.0868, 0.0868, 0.1081, 0.1081};

    // Precipitation lapse rate derived from stations located in the Upper Bow River basin by Logan Fang :
    //   Stations Banff VS Sunshine; Station Lake Louise VS Skoki
    // (used in CRHM for Marmot Creek domain)
    //(100 m)^-1
   // double monthly_factors[] = {0.1932,0.1932 ,0.1932, 0.1182, 0.1182, 0.1182, 0.0592, 0.0592, 0.09292, 0.09292, 0.1932, 0.1932};


    for(auto& mf: monthly_factors)
    {
        mf /= 100.0; //to m^-1
    }
    std::vector< boost::tuple<double, double, double> > ppt;
    std::vector< boost::tuple<double, double, double> > staion_z;
    for (auto& s : face->stations())
    {
        if( is_nan((*s)["p"_s]))
            continue;
        double p = (*s)["p"_s];
        ppt.push_back( boost::make_tuple(s->x(), s->y(), p ) );
        staion_z.push_back( boost::make_tuple(s->x(), s->y(), s->z() ) );
    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double p0 = face->get_module_data<data>(ID)->interp(ppt, query);
    double z0 = face->get_module_data<data>(ID)->interp(staion_z,query);
    double z = face->get_z();
    double slp = face->slope();

    int month = global_param->month()-1;
    double lapse = monthly_factors[month];
    double adj_fac = 1. + lapse*(z-z0);

    // Limit the precipitation-elevation adjustment factor in the range 0.5 - 1.5
    adj_fac = std::max(0.5,adj_fac);
    adj_fac = std::min(1.5,adj_fac);
    
    double P = p0*adj_fac;
    
    double P_fin;

   // Correct precipitation input using triangle slope when input preciptation are given for the horizontally projected area.
    if(apply_cosine_correction)
    {
        P_fin = Atmosphere::corr_precip_slope(P,slp);
    } else
    { 
        P_fin = P;
    }

    P_fin =  std::max(0.0,P_fin);
    (*face)["p"_s]= P_fin;

    P = std::max(0.0,P);
    (*face)["p_no_slope"_s]= P;

}

p_lapse::~p_lapse()
{

}
