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

//
// Created by Donovan Allum 2025
//

#pragma once

#include "global.hpp"
#include "triangulation.hpp"
#include <limits>
#include <optional>
#include <boost/shared_ptr.hpp>
#include <boost/property_tree/ptree.hpp>
#include <cstddef>
#include <stdexcept>
#include <cassert>
#include <cstdint>
#include <concepts>
#include <type_traits>

/**
 * @brief Base class for data management with lazy caching and output handling.
 *
 * This class provides a framework for managing data with optional caching,
 * lazy evaluation, and output assignment. It uses CRTP-like patterns with
 * concepts to enforce interface contracts at compile-time.
 *
 * Key features:
 * - Automatic cache initialization and staleness checking
 * - Lazy evaluation of values with cache-aware updates
 * - Type-safe output assignment with runtime checks
 * - Compile-time concept enforcement for derived cache types and value/output providers
 *
 * The class is designed for computational scenarios where:
 * - Data may be expensive to compute and should be cached
 * - Cache validity depends on external timestep counters
 * - Values need to be fetched lazily when accessed
 * - Outputs must be accumulated safely
 *
 * Usage example:
 * @code
 * struct MyCache : public cache_base {
 *     double temperature = default_value<double>();
 *     int iteration_count = default_value<int>();
 * };
 *
 * class MyData : public data_base<MyCache> {
 * public:
 *     MyData(const mesh_elem& face_in, const boost::shared_ptr<global> param, 
 *            const pt::ptree& cfg) : data_base(face_in, param, cfg) {}
 *     
 *     void compute_temperature() {
 *         update_value(
 *             [this]() -> auto& { return cache_->temperature; },
 *             [this]() { return expensive_temperature_calculation(); }
 *         );
 *     }
 *     
 *     void add_to_output(double contribution) {
 *         set_output(
 *             [this]() -> auto& { return output_variable; },
 *             contribution
 *         );
 *     }
 * };
 * @endcode
 *
 * @tparam CacheType A type derived from cache_base that provides storage
 *         for cached values. Must satisfy the CacheRules concept.
 */

struct cache_base 
{
    int64_t last_timestep = -1;

    template<typename T>
    static constexpr T default_value() {
        if constexpr (std::is_floating_point_v<T>)
            return std::numeric_limits<T>::quiet_NaN();
        else if constexpr (std::is_integral_v<T>)
            return std::numeric_limits<T>::min();
        else 
            return T{};
    };

};

namespace data_base_concepts {
    template<typename C>
    concept CacheRules = std::derived_from<C,cache_base>;

    template<typename V>
    concept ValueRules = 
        requires(V v) {
        {v()} -> std::same_as<std::add_lvalue_reference_t<decltype(v())>>;
    };

    template<typename O,typename T>
    concept OutputRules = ValueRules<O> &&
    requires(O o,const T t)
    {
        {o()} -> std::convertible_to<T>;
        {o() += t};
    };
};
    
namespace pt = boost::property_tree;

template<data_base_concepts::CacheRules CacheType>
class data_base {
    
    template<typename T>
    bool constexpr is_unset(T t)
    {
        if constexpr (!std::is_floating_point_v<T>)
            return t == cache_base::default_value<T>();
        else 
            return std::isnan(t);
    };

    bool is_stale();

    void init_cache();

protected:
    
    data_base(const mesh_elem& face_in, const boost::shared_ptr<global> param, 
            const pt::ptree& cfg, const bool istest = false);
    ~data_base() {};
    
    const mesh_elem face{nullptr};
    const boost::shared_ptr<global> global_param;
    const pt::ptree& cfg_;
    mutable std::optional<CacheType> cache_;

    template<data_base_concepts::ValueRules Value,typename Fetch>
    void update_value(Value&& value, const Fetch& fetch);

    template<typename T,data_base_concepts::OutputRules<T> Output>
    void set_output(Output&& output,const T t);

public:
    void reset_cache() { cache_.reset(); };
    const std::optional<CacheType>& get_cache() { return cache_; }; 
};

template<data_base_concepts::CacheRules CacheType>
void data_base<CacheType>::init_cache() {
    if (!cache_ || is_stale()) {
        cache_.emplace();
        cache_->last_timestep = global_param->timestep_counter;
    }
}

template<data_base_concepts::CacheRules CacheType>
bool data_base<CacheType>::is_stale()
{
    return cache_->last_timestep != global_param->timestep_counter;
}

template<data_base_concepts::CacheRules CacheType>
data_base<CacheType>::data_base(const mesh_elem& face_in, const boost::shared_ptr<global> param, 
        const pt::ptree& cfg, const bool istest) : face(face_in), global_param(param), cfg_(cfg)
{
	// Optional istest parameter only exists to skip these tests during tests of this class where we aren't testing whether the face object has been set correctly.
	// This means tests show that the underlying functions work as intended
    if (!face->is_valid() && !istest)
        throw std::invalid_argument("Face handle points to an invalid face");

    if (!global_param && !istest)
        throw std::invalid_argument("global parameter holder is null");
};

template<data_base_concepts::CacheRules CacheType>
template<data_base_concepts::ValueRules Value,typename Fetch>
void data_base<CacheType>::update_value(Value&& value, const Fetch& fetch) {
    init_cache();
    
    auto& V = value(); 
    if ( is_unset(V) )
    {
        V = fetch();
    }

};

template<data_base_concepts::CacheRules CacheType>
template<typename T,data_base_concepts::OutputRules<T> Output>
void data_base<CacheType>::set_output(Output&& output,const T t)
{
    init_cache();

    output() += t;
};
