// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>

namespace snowslide
{
enum class Profile
{
    bernhardt_schulz,
    chm_vionnet
};

// Metres of vertical snow depth digitized every 5 degrees from Bernhardt & Schulz
// (2010), Fig. 2b. The printed curve has about 0.17 m per pixel on its depth axis.
// Section 3 specifies the >=75 degree plateau separately, in SWE.
inline double paper_depth(double degrees)
{
    constexpr std::array depths{5.0, 3.4, 2.4, 1.8, 1.4, 1.1, 0.9, 0.8, 0.6, 0.6, 0.5};
    const double index = std::clamp((degrees - 25.0) / 5.0, 0.0, 10.0);
    const size_t i = std::min<size_t>(static_cast<size_t>(index), depths.size() - 2);
    return std::lerp(depths[i], depths[i + 1], index - i);
}

inline double chm_depth(double radians, double canopy, double coefficient, double exponent)
{
    const double degrees = radians * 180.0 / std::numbers::pi;
    return std::min(
        20.0, std::max(coefficient * std::pow(std::max(10.0, degrees), exponent) * std::max(0.01, std::cos(radians)),
                       canopy));
}

// SWE and signed drift must use the same units. Wind deposits enter the snow column after routing.
inline void reserve_wind_erosion(double drift, double& swe, double& normal, double& vertical)
{
    if (drift < 0 && swe > 0)
    {
        const double remaining = std::max(0.0, swe + drift);
        const double fraction = remaining / swe;
        normal *= fraction;
        vertical *= fraction;
        swe = remaining;
    }
}

struct Release
{
    double mass;           // kg/m2
    double normal_removed; // m
    double vertical_left;  // m
    double capacity_left;  // kg/m2 of new deposits that can stay on the track
    double max_depth;      // equivalent slope-normal threshold (m)
};

inline Release release(Profile profile, double radians, double vegetation, double coefficient, double exponent,
                       double mass, double normal, double vertical, double deposit_density)
{
    const double cosine = std::max(0.001, std::cos(radians));
    const double degrees = radians * 180.0 / std::numbers::pi;
    const bool paper = profile == Profile::bernhardt_schulz;
    // Avoid converting a steep-slope SWE floor to depth with a fixed snow density.
    const double normal_density = normal > 0 ? mass / normal : deposit_density;
    const double vertical_density = vertical > 0 ? mass / vertical : deposit_density * cosine;
    const double normal_cap = paper ? 0.0 : chm_depth(radians, vegetation, coefficient, exponent);
    const double vertical_cap = paper ? std::max(vegetation, degrees < 75.0 ? paper_depth(degrees) : 0.0) : 0.0;
    const double capacity = paper ? std::max(50.0, vertical_cap * vertical_density) : normal_cap * normal_density;
    const bool eligible = !paper || radians > 25.0 * std::numbers::pi / 180.0;
    const double moved = eligible ? std::max(0.0, mass - capacity) : 0.0;
    const double fraction = mass > 0 ? moved / mass : 0.0;
    const double vertical_left = vertical * (1.0 - fraction);
    double room = (normal_cap - normal * (1.0 - fraction)) * deposit_density;
    if (paper)
    {
        room = std::max(50.0 - (mass - moved), (vertical_cap - vertical_left) * deposit_density * cosine);
    }
    return {moved, normal * fraction, vertical_left,
            eligible ? std::max(0.0, room) : std::numeric_limits<double>::infinity(),
            paper ? capacity / normal_density : normal_cap};
}

struct Flux
{
    double value, donor_derivative, neighbour_derivative;
};

// Positive out of i. The Jacobian holds the nonorthogonal correction fixed.
inline Flux edge_flux(double mi, double mj, double head_drop, double transmissibility, double mobile_density)
{
    const double velocity = transmissibility * head_drop;
    const double upwind = velocity >= 0.0 ? mi : mj;
    const double pressure = transmissibility * upwind / mobile_density;
    return {velocity * upwind, std::max(0.0, velocity) + pressure, std::min(0.0, velocity) - pressure};
}

inline double deposit(double mobile, double step, double room, bool has_outflow)
{
    const double settling = -std::expm1(-step) * mobile;
    return has_outflow ? std::min(settling, std::max(0.0, room)) : settling;
}
} // namespace snowslide
