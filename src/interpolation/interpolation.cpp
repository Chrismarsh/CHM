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

#include "interpolation.h"

interpolation::interpolation(interp_alg ia, size_t size,
                             std::map<std::string,std::string> config)
{
    init(ia,size,config);
}

void interpolation::init(interp_alg ia, size_t size,std::map<std::string,std::string> config)
{

    if (ia == interp_alg::tpspline)
    {
        base = boost::make_shared<thin_plate_spline>(size,config);
    }
    else if(ia == interp_alg::idw)
    {
        base = boost::make_shared<inv_dist>();
    }
    else if(ia == interp_alg::nearest_sta)
    {
        base = boost::make_shared<nearest>();
    }
    else
    {
        BOOST_THROW_EXCEPTION(interp_unknown_type() << errstr_info("Unknown interpolation type"));
    }


    this->size = size;
    this->ia = ia;
}

interpolation::interpolation()
{
    base = nullptr;
}
interpolation::~interpolation()
{

}

double interpolation::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)
{
    if (sample_points.size() == 0)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info("Interpolation sample point length = 0."));
    }

    if (sample_points.size() > 15 && ia == interp_alg::tpspline)
    {
        LOG_WARNING << "More than 15 sample points is likely to cause slow downs";
    }

//    if(sample_points.size() != this->size || this->size == 0)
//    {
//        this->init(ia,sample_points.size());
//    }

    return base->operator()(sample_points,query_point);
}
