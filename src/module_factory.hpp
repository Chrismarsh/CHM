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

#pragma once

#include "module_base.hpp"
#include "logger.hpp"
#include "exception.hpp"

#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
namespace pt = boost::property_tree;

// Modules to build against
#include "iswr.hpp"
#include "Marsh_shading_iswr.hpp"
#include "const_llra_ta.hpp"
#include "Kunkel_monthlyTd_rh.hpp"
#include "Liston_monthly_llra_ta.hpp"
#include "Sicart_ilwr.hpp"
#include "Liston_wind.hpp"
#include "PenmanMonteith_evaporation.hpp"
#include "Thornton_p.hpp"
#include "Walcek_cloud.hpp"
#include "Harder_precip_phase.hpp"
#include "Burridge_iswr.h"
#include "Iqbal_iswr.h"
#include "iswr_from_obs.h"
#include "Dodson_NSA_ta.h"
#include "Thornton_p.hpp"
#include "Thornton_var_p.h"
#include "rh_from_obs.h"
#include "kunkel_rh.hpp"
#include "snobal.hpp"
#include "Richard_albedo.hpp"
#include "snowpack.hpp"
#include "point_mode.hpp"
#include "threshold_p_phase.hpp"
#include "Gray_inf.hpp"
#include "Longwave_from_obs.hpp"
#include "Simple_Canopy.hpp"
#include "Dist_tlapse.hpp"
#include "scale_wind_vert.hpp"
#include "solar.hpp"
#include "rh_no_lapse.hpp"
#include "t_no_lapse.hpp"
#include "p_no_lapse.hpp"
#include "lw_no_lapse.hpp"
#include "uniform_wind.hpp"
#include "fast_shadow.hpp"
#include "deform_mesh.hpp"
#include "crop_rotation.hpp"
#include "sub_grid.hpp"
#include "t_monthly_lapse.hpp"
#include "snow_slide.hpp"
#include "PBSM3D.hpp"
#include "fetchr.hpp"
#include "MS_wind.hpp"
#include "Cullen_monthly_llra_ta.hpp"
class module_factory
{
public:
    module_base* get(std::string ID, pt::ptree config);
};
