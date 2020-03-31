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


#include "iswr_from_nwp.hpp"
REGISTER_MODULE_CPP(iswr_from_nwp);


iswr_from_nwp::iswr_from_nwp(config_file cfg)
        : module_base("iswr_from_nwp", parallel::data, cfg)
{
    depends_from_met("Qsi");
    depends_from_met("Qsi_diff");

    provides("iswr_direct_no_slope");
    provides("iswr_diffuse_no_slope");
    provides("iswr_observed");

}
iswr_from_nwp::~iswr_from_nwp()
{

}
void iswr_from_nwp::init(mesh& domain)
{

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {

	       auto face = domain->face(i);
	       auto d = face->make_module_data<data>(ID);
	       d->interp.init(global_param->interp_algorithm,face->stations().size() );

    }

}
void iswr_from_nwp::run(mesh_elem &face)
{
    //interpolate all the measured qsi and qsi_diff from the NWP model

    //lower all the station values to sea level prior to the interpolation
    std::vector< boost::tuple<double, double, double> > lowered_values;
    std::vector< boost::tuple<double, double, double> > lowered_values2;
    for (auto& s : face->stations())
    {
        if( (is_nan((*s)["Qsi"_s])) || (is_nan((*s)["Qsi_diff"_s])))
            continue;
        double v = (*s)["Qsi"_s];
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), v ) );
        double vv = (*s)["Qsi_diff"_s];
        lowered_values2.push_back( boost::make_tuple(s->x(), s->y(), vv ) );
    }

    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    // Read interpolated total and diffuse iswr
    double iswr_observed =face->get_module_data<data>(ID)->interp(lowered_values, query);
    double split_diff =face->get_module_data<data>(ID)->interp(lowered_values2, query);

    // Compute direct part
    double split_dir = iswr_observed - split_diff;

    split_dir = std::max(0.0,split_dir);
    split_diff = std::max(0.0,split_diff);

    (*face)["iswr_direct_no_slope"_s]=split_dir;
    (*face)["iswr_diffuse_no_slope"_s]=split_diff;
    (*face)["iswr_observed"_s]=iswr_observed;

}
