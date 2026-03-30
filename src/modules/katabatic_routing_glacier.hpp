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

#include "triangulation.hpp"
#include "module_base.hpp"
#include "Glacier.hpp"
#include "katabatic_melt_energy.hpp"
#include "melt_routing_glacier.hpp"
#include "data_base.hpp"
#include "PhysConst.h"

/**
 * \ingroup modules infil soils exp
 * @{
 * \class Gray_inf
 *
 *
 * Estimates areal snowmelt infiltration into frozen soils for:
 *    a) Restricted -  Water entry impeded by surface conditions
 *    b) Limited - Capiliary flow dominates and water flow influenced by soil physical properties
 *    c) Unlimited - Gravity flow dominates
 *
 * **Depends:**
 * - Snow water equivalent "swe" [mm]
 * - Snow melt for interval "snowmelt_int" [\f$mm \cdot dt^{-1}\f$]
 *
 * **Provides:**
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
 * .. note::
 *    Has hardcoded soil parameters that need to be read from the mesh parameters.
 *
 * \endrst
 *
 * **References:**
 * - Gray, D., Toth, B., Zhao, L., Pomeroy, J., Granger, R. (2001). Estimating areal snowmelt infiltration into frozen soils
 * Hydrological Processes  15(16), 3095-3111. https://dx.doi.org/10.1002/hyp.320
 * @}
 */
class katabatic_routing_glacier : public module_base
{
REGISTER_MODULE_HPP(katabatic_routing_glacier)
public:
    katabatic_routing_glacier(config_file cfg);

    ~katabatic_routing_glacier();

    void run(mesh_elem &face) override;
    void init(mesh& domain) override;

	struct Cache : public cache_base
	{
		Cache() {};
		// All submodule inputs
		double swe = cache_base::default_value<double>();
        double melt_energy = cache_base::default_value<double>();
		double glacier_temperature = cache_base::default_value<double>();
        double update_now = cache_base::default_value<double>();
		double air_temperature = cache_base::default_value<double>();
		double air_pressure = cache_base::default_value<double>();
		double vapour_pressure = cache_base::default_value<double>();
		double vapour_pressure_surface = cache_base::default_value<double>();
		double lapse_rate = cache_base::default_value<double>();
		double snowmelt = cache_base::default_value<double>();

		// All submodule outputs
    	double glacier_water_equivalent = 0.0;
        double total_depth = 0.0;
        double firnmelt = 0.0;
        double icemelt = 0.0;
		double latent_heat = 0.0;
		double sensible_heat = 0.0;
		double snowmelt_delayed = 0.0;
		double firnmelt_delayed = 0.0;
		double icemelt_delayed = 0.0;
		double total_delayed = 0.0;
	};

    class data : public face_info, public data_base<Cache>
    {
    public:
		data(const mesh_elem& face, const std::shared_ptr<global> param, const config_file cfg,
				const Glacier::Params* g_p, const GlacierRouting::Params* r_p); 
		
		GlacierRouting::State routing_state;
		Glacier::State glacier_state;

		Units::Kelvin glacier_temperature();
		Units::Kelvin air_temperature();
		Units::Pa air_pressure();
		Units::Pa vapour_pressure();
		Units::Pa vapour_pressure_surface();
		Units::LapseRateSI lapse_rate();
		void latent_heat(const double v);
		void sensible_heat(const double v);

        const Units::Milimetres snowmelt();
        const Units::Milimetres firnmelt();
        const Units::Milimetres icemelt();
        void snowmelt_delayed(double v);
        void firnmelt_delayed(double v);
        void icemelt_delayed(double v);
		void total_delayed(double v);

		const Units::Milimetres swe();
		const Units::Watts_per_m2 melt_energy();
		bool update_now();

		void glacier_water_equivalent(const double out);
		void total_depth(const double out);
		void firn_melt(const double out);
		void ice_melt(const double out);

		double firn_emissivity = 0.0;
		double ice_emissivity = 0.0;
		double total_energy = 0.0;
	};

	// Adapter that satisfies KatabaticData concept
	class katabatic_view
	{
		data& d;
		mesh_elem& face;
	public:
		explicit katabatic_view(data& d,mesh_elem& face)
			: d(d), face(face){};
		~katabatic_view();
		Units::Kelvin glacier_temperature() { return d.glacier_temperature(); }
		Units::Kelvin air_temperature() { return d.air_temperature(); }
		Units::Pa air_pressure() { return d.air_pressure(); }
		Units::Pa vapour_pressure() { return d.vapour_pressure(); }
		Units::Pa vapour_pressure_surface() { return d.vapour_pressure_surface(); }
		Units::LapseRateSI lapse_rate() { return d.lapse_rate(); }
		void latent_heat(const double v) { d.latent_heat(v); }
		void sensible_heat(const double v) { d.sensible_heat(v); }
	};


	// Adapter that satisfies GlacierRoutingData concept
    class routing_view {
        data& d;
		mesh_elem& face;
    public:
        explicit routing_view(data& d,mesh_elem& face) : d(d),face(face) {}
		~routing_view();
		GlacierRouting::State& get_state() { return d.routing_state; }  // no name collision!
        const Units::Milimetres snowmelt() { return d.snowmelt(); }
        const Units::Milimetres firnmelt() { return d.firnmelt(); }
        const Units::Milimetres icemelt() { return d.icemelt(); }
        void snowmelt_delayed(const double v) { d.snowmelt_delayed(v); }
        void firnmelt_delayed(const double v) { d.firnmelt_delayed(v); }
        void icemelt_delayed(const double v) { d.icemelt_delayed(v); }
        void total_delayed(const double v) { d.total_delayed(v); }
    };

	// Similarly for Glacier::Model's concept
    class glacier_view {
        data& d;
		mesh_elem& face;
    public:
        explicit glacier_view(data& d,mesh_elem& face) 
			: d(d), face(face) {}
		~glacier_view();
		Glacier::State& get_state() { return d.glacier_state; }
		const Units::Milimetres swe() { return d.swe(); }
		const Units::Watts_per_m2 melt_energy() { return Units::Watts_per_m2{d.total_energy}; }
		const Units::Celsius glacier_temperature() { return d.glacier_temperature(); }
		bool update_now() { return d.update_now(); }

		void glacier_water_equivalent(const double out) { d.glacier_water_equivalent(out); }
		void total_depth(const double out) { d.total_depth(out); }
		void firn_melt(const double out) { d.firn_melt(out); }
		void ice_melt(const double out) { d.ice_melt(out); }
    };

private:
    katabatic_melt_energy::Model<katabatic_view> katabatic;
    GlacierRouting::Model<routing_view> routing;
    Glacier::Model<glacier_view> glacier;

	bool is_new_day();
	double rain_sun_energy(const mesh_elem& face);
	double rain_sun_energy(const mesh_elem& face,data& d);
};
