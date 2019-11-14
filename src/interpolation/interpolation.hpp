/* * Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
 * modular unstructured mesh based approach for hydrological modelling
 * Copyright (C) 2018 Christopher Marsh
 *
 * This file is part of Canadian Hydrological Model.
 *
 * Canadian Hydrological Model is free software: you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Canadian Hydrological Model is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Canadian Hydrological Model.  If not, see
 * <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "interp_base.hpp"
#include "inv_dist.hpp"
#include "nearest.hpp"
#include "TPSpline.hpp"

#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "logger.hpp"
enum interp_alg
{
    tpspline,
    idw,
    nearest_sta
};

class interpolation
{
public:
    /*
     * Some of the interpolators need to allocate memory, and if called multiple times, it is faster to preinit and set size
     * as required where size is the number of items used to interpolate from. E.g., # stations.
     */
    interpolation(interp_alg ia, size_t size=0,
                  std::map<std::string,std::string> config = std::map<std::string,std::string>());
    interpolation();

    ~interpolation();

    void init(interp_alg ia, size_t size=0, std::map<std::string,std::string> config = std::map<std::string,std::string>());

    double operator()(std::vector< boost::tuple<double,double,double> >& sample_points, boost::tuple<double,double,double>& query_point);
    boost::shared_ptr<interp_base> base;
private:

    size_t size;
    interp_alg ia;
};
