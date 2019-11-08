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

// hash functions
#include "utility/BBhash.h"
#include "utility/wyhash.h"

#include <string>
#include <vector>

class variablestorage
{
  public:
    variablestorage();
    variablestorage(std::set<std::string>& variables);
    ~variablestorage();

    double& operator[](const uint64_t& variable);
    double& operator[](const std::string& variable);

    bool has(const uint64_t& hash);
    bool has(const std::string& variable);

    std::vector<std::string> variables();

    void init(std::set<std::string>& variables);

  private:
    const uint64_t seed = 2654435761U;
    template <typename Item> class wyandFunctor
    {
      public:
        uint64_t operator ()  (const Item& key, uint64_t seed = seed) const
        {
            return wyhash(&key, sizeof(Item), seed);
        }

    };
    typedef wyandFunctor<uint64_t> hasher_t;
    typedef boomphf::mphf< uint64_t, hasher_t  > boophf_t;

    // Note that we have to explicitly check if what we get back is what we wanted as
    // mphf do not guarantee what asking for something outside of the map returns a sane answer
    // https://github.com/rizkg/BBHash/issues/12

    // perfect hashfn + variable storage
    std::unique_ptr<boophf_t> _variable_bphf;
    std::vector<var> _variables;

    size_t _size;

};

