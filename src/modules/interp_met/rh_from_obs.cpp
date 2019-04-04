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


#include "rh_from_obs.hpp"
REGISTER_MODULE_CPP(rh_from_obs);

rh_from_obs::rh_from_obs(config_file cfg)
  : module_base("rh_from_obs", parallel::data, cfg)
{
    depends_from_met("rh");
    depends_from_met("t");
    depends("t");
    provides("rh");
}
rh_from_obs::~rh_from_obs()
{

}
void rh_from_obs::init(mesh& domain)
{
    ompException oe;
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
      oe.Run([&]
	     {
	       auto face = domain->face(i);
	       auto d = face->make_module_data<data>(ID);
	       d->interp.init(global_param->interp_algorithm,global_param->get_stations( face->get_x(), face->get_y()).size());
	     });
    }
    oe.Rethrow();
}

void rh_from_obs::run(mesh_elem& face)
{
    //generate lapse rates
    std::vector<double> sea;
    std::vector<double> sz;

    static boost::posix_time::ptime last_update;
    static double lapse=-999.0;

    //we have to run this at least once, use the first bool flag to do so
    //otherwise just check that the last update was at a different time step, if so calculate the lapse rate.
    //otherwise, just used the stored lapse rate
    if(last_update != global_param->posix_time() )
    {
        for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
        {
            if( is_nan(s->get("t")) || is_nan(s->get("rh")))
                continue;
            double rh = s->get("rh")/100.;
            double t = s->get("t");
            double es = mio::Atmosphere::vaporSaturationPressure(t+273.15);
            double ea = rh * es;
            sea.push_back( ea  );
            sz.push_back( s->z());
        }


        // least squares linear fit to these points ( p v. z)
        double c0, c1, cov00, cov01, cov11, chisq;

        gsl_fit_linear (&sz[0], 1,&sea[0], 1, sz.size(),
                        &c0, &c1, &cov00, &cov01, &cov11,
                        &chisq);
        lapse = c1;//use the slope (y=mx+b c1 == m)
        last_update = global_param->posix_time();
    }

    std::vector< boost::tuple<double, double, double> > lowered_values;
    for (auto& s : global_param->get_stations( face->get_x(), face->get_y()))
    {
        if( is_nan(s->get("t")) || is_nan(s->get("rh")))
            continue;

        double rh = s->get("rh")/100.;
        double t = s->get("t");
        double es = mio::Atmosphere::vaporSaturationPressure(t+273.15);
        double ea = rh * es;
        double z = s->z();
        ea = ea + lapse*(0.0-z);
        lowered_values.push_back( boost::make_tuple(s->x(), s->y(), ea ) );

    }


    auto query = boost::make_tuple(face->get_x(), face->get_y(), face->get_z());
    double ea = face->get_module_data<data>(ID)->interp(lowered_values, query);

    //raise it back up
    ea = ea + lapse*( face->get_z() - 0.0);

    double es = mio::Atmosphere::vaporSaturationPressure((*face)["t"_s]+273.15);
    double rh = ea/es*100.0;

    rh = std::min(rh,100.0);
    rh = std::max(10.0,rh);

    (*face)["rh"_s]=rh;


}
