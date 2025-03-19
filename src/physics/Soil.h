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

//
// Created by nicway on 12/07/16.
//

#pragma once

#include <iostream>
#include <sparsehash/dense_hash_map>
#include "exception.hpp"

namespace Soil {
    /********* Soil ************/

    enum GATable {PSI, KSAT, WILT, FCAP, PORG, PORE, AIRENT, PORESZ, AVAIL}; // Used for mapping the soil table, PSI and KSAT are used, the others are unused but may but used in the future or other modules.    
    // see GreenAmpt module in CRHM wiki for details.
    //

    using mymap = google::dense_hash_map<std::string, double>;

    class _soils_base
    {
    public:
        virtual double porosity(std::string) const = 0;
        virtual double pore_size_dist(std::string) const = 0;
        virtual double wilt_point(std::string) const = 0;
        virtual double air_entry_tension(std::string) const = 0;
        virtual double capillary_suction(std::string) const = 0;
        virtual double saturated_conductivity(std::string) const = 0;

        virtual double ayers_texture(std::string texture, std::string ground_cover) const = 0;
        double lookup(const mymap& map, const std::string key) const;
        virtual ~_soils_base() = default;

    private:
        virtual void _make_hash() = 0;
    };

    class soils_na : public _soils_base
    {
    public:
        soils_na();
        ~soils_na();

        double porosity(std::string soil_type) const override;
        double pore_size_dist(std::string soil_type) const override;
        double wilt_point(std::string soil_type) const override;
        double air_entry_tension(std::string soil_type) const override;
        double capillary_suction(std::string soil_type) const override;
        double saturated_conductivity(std::string soil_type) const override;
        double ayers_texture(std::string texture, std::string ground_cover) const override;


    private:

        mymap _porosity;
        mymap _pore_size_dist;
        mymap _wilt_point;
        mymap _air_entry_tension;
        mymap _capillary_suction;
        mymap _saturated_conductivity;

        google::dense_hash_map< std::string, mymap> _ayers_texture;
        
        void _make_hash() override;


    };

};
