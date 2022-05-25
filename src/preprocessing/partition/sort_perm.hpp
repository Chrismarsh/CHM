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
#include <tbb/parallel_sort.h>

// adapted from
// https://stackoverflow.com/a/17074810/410074
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    Compare compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
              [&](std::size_t i, std::size_t j){ return compare(vec.at(i), vec.at(j)); });
    return p;
}

template <typename T>
void apply_permutation_in_place(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<bool> done(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        if (done.at(i))
        {
            continue;
        }
        done.at(i) = true;
        std::size_t prev_j = i;
        std::size_t j = p.at(i);
        while (i != j)
        {
            std::swap(vec.at(prev_j), vec.at(j));
            done.at(j) = true;
            prev_j = j;
            j = p.at(j);
        }
    }
}

template <typename T>
std::vector<T> apply_permutation_transform(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
                   [&](std::size_t i){ return vec[i]; });

    return sorted_vec;
}

template <typename T>
void apply_permutation_in_place_naive(
    std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    assert(vec.size() == p.size());
    std::vector<T> res(vec.size());
    for (std::size_t i = 0; i < vec.size(); ++i)
    {
        res.at(p.at(i)) = vec.at(i);
    }
    std::swap(vec,res);
}