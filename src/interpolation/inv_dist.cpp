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

#include "inv_dist.hpp"

inv_dist::inv_dist()
{
    
}

inv_dist::~inv_dist()
{
    
}

double inv_dist::operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point)
{
    
    double numerator = 0.0;
    double denominator = 0.0;

    double z0 = 0;
    if (sample_points.size() == 0)
    {
        BOOST_THROW_EXCEPTION( interpolation_error()
                                << errstr_info("IDW requires >=1 stations"));
    }
    for(size_t i=0;i<sample_points.size();i++)
    {
        double z = sample_points.at(i).get<2>();

        double sx = sample_points.at(i).get<0>();
        double sy = sample_points.at(i).get<1>();
        
        double ex = query_point.get<0>();
        double ey = query_point.get<1>();
        
        double xdiff = (sx  - ex);
        double ydiff = (sy  - ey);
        double di = pow(sqrt(
                             pow(xdiff,2.0) + pow(ydiff,2.0)
                            )
                        ,2.0
                    );
        if(di == 0)
        {
                numerator = z;
                denominator = 1.0;
        }
        else
        {
                numerator += z/ di;
                denominator += 1 / di; 
        }
    }

   z0 = (numerator/denominator);

   return z0;

}
