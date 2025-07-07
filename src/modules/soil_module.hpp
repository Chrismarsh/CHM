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

#include "logger.hpp"
#include "triangulation.hpp"
#include "module_base.hpp"
#include "TPSpline.hpp"
#include <cmath>
#include <memory>
#include "Soil.h"
#include "soil_two_layer.hpp"
#include "soil_ET.hpp"
#include "K_estimate.hpp"
#include "XG_algorithm.hpp"

/**
 * \ingroup modules infil soil_module exp
 * @{
 * \class soil_module *
 * 
 * Organizes soil process models (separate classes)
 *      - estimating freeze/thaw fronts using the XG algorithm: XG-freeze_thaw
 *      - layer per-time-step outflow: k_estimate
 *      - two layer soil model: soil_two_layer
 *      - soil Evapotranspiration: soil_ET
 *
 *      Each are found in soil_submodules/
 *
 * **Depends:**
 * - Snow water equivalent "swe" [mm]
 * - thaw front in soil "thaw_front_depth" [mm]
 * - freeze front in soil "freeze_front_depth" [mm]
 * - potential evapotranspiration "potential_ET" [mm]
 * - infiltrated moisture "inf" [mm]
 * - runoff/excess from infiltration scheme "runoff" [mm]
 * - leftover from routing between triangles "routing_residual" [mm]             
 *
 * **Provides:** TODO match provides with actual .cpp file
 * - Infiltration "inf" [\f$mm \cdot dt^{-1}\f$]
 * - Total infiltration "total_inf" [mm]
 * - Total infiltration excess "total_excess" [mm]
 * - Total runoff "runoff" [mm]
 * - Total soil storage "soil_storage"
 * - Potential infiltration "potential_inf"
 * - Opportunity time for infiltration to occur "opportunity_time"
 * - Available storage for water of the soil "available_storage"
 *
 * \rst
 * .. note:: TODO add any notes
 *    Has hardcoded soil parameters that need to be read from the mesh parameters.
 *
 * \endrst
 *
 * **References:** TODO add any references
 * - Gray, D., Toth, B., Zhao, L., Pomeroy, J., Granger, R. (2001). Estimating areal snowmelt infiltration into frozen soils
 * Hydrological Processes  15(16), 3095-3111. https://dx.doi.org/10.1002/hyp.320
 * @}
 */

class soil_module : public module_base
{
REGISTER_MODULE_HPP(soil_module);
public:
    soil_module(config_file cfg);

    ~soil_module();

    void run(mesh_elem& face);
    void init(mesh& domain);

    class data : public face_info, public two_layer_DTO, public soil_ET_DTO
    {
    public:
        std::unique_ptr<soil_base> soil_layers;
        std::unique_ptr<soil_base> ET;
        // custom deletor that does nothing to make sure it doesn't try to delete the mesh_elem it points to
        // REMOVED same reason as is_lake below, may 2025
		//mesh_elem* my_face;//(nullptr, [](mesh_elem* ptr) {});
        std::unique_ptr<I_K_estimate> K_estimator;
        // overridden
        // REMOVED is_lake as a function and is now a variable storing the result of is_water (May 2025)
		//bool is_lake(soil_ET_DTO& DTO) override;
        int get_dt() override;
        bool get_new_day() override;
        // custom deletor that does nothing to make sure it doesn't try to delete the mesh_elem it points to
        soil_module* local_module;//(nullptr, [](soil_module*) {});
    
        std::unique_ptr<XG_algorithm::params> P;
        std::unique_ptr<XG_algorithm::state> S;

        bool first_day = true;
    };

    void set_local_module(soil_module::data& d);


private:
    
    std::unique_ptr<Soil::_soils_base> SoilDataObj;

    void get_soil_inputs(mesh_elem& face, data& d);
    void set_soil_outputs(mesh_elem& face, data& d);
    void set_soil_params(mesh_elem& face, soil_module::data& d);
    void set_ET_params(mesh_elem& face, soil_module::data& d);
    void initial_soil_conditions(mesh_elem& face, soil_module::data& d);

    int compare_substring(std::string& type, std::string sub);

    XG_algorithm get_XG(mesh_elem& face,data& d);
    void init_param_state_XG(mesh_elem& face,data& d);
    bool is_new_day();
    struct XG_shared_const
    {
        double Trigthrhld;
        size_t num_layers;
        double theta_min;
        size_t freeze_kw_ki_update;
        size_t thaw_ki_kw_update;
        size_t k_update;
        size_t time_step_per_day;
        bool calc_conductivity;
    };
    XG_shared_const C;
};
