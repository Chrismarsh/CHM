// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include "crocus_bridge.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace crocus
{
// Remove surface material before the vertical step. Enthalpy is extensive; the
// density, grain properties and age of a partially eroded layer are intensive.
inline bool erode(int n, float* state, float& trace, double demand)
{
    if (demand <= 0)
        return true;
    int active = 0;
    double available = trace;
    while (active < n && state[active] > 0)
        available += state[active++];
    const double tolerance = 4 * std::numeric_limits<float>::epsilon() * std::max(1., available);
    if (demand > available + tolerance)
        return false;

    const double from_trace = std::min(demand, static_cast<double>(trace));
    trace -= from_trace;
    demand -= from_trace;
    int first = 0;
    while (first < active && demand >= state[first])
        demand -= state[first++];
    if (first < active && demand > 0)
    {
        const float old = state[first];
        state[first] = static_cast<float>(old - demand);
        state[heat * n + first] *= state[first] / old;
    }
    if (first)
    {
        for (int f = 0; f < field_count; ++f)
        {
            auto* row = state + f * n;
            std::move(row + first, row + active, row);
            std::fill(row + active - first, row + active, 0.f);
        }
    }
    return true;
}

inline double avalanche_mass_per_area(double water_volume, double area) { return 1000. * water_volume / area; }
} // namespace crocus
