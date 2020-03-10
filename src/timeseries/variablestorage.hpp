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
#include "utility/xxh64.hpp"


#include "logger.hpp"
#include "exception.hpp"

#include <string>
#include <vector>
#include <set>
#include <type_traits>


template<typename T = double>
class variablestorage
{
  public:
    variablestorage();

    /// Initialize the storage with a set of variables. Values default to -9999
    /// @param variables
    variablestorage(std::set<std::string>& variables);
    ~variablestorage();

    /// Get and set the variable to a specific value. Use _s for compile-time hash.
    /// Throws if not found or init/ctor not yet called.
    /// @param variable
    /// @return
    T& operator[](const uint64_t& variable);
    /// Get and set the variable to a specific value.
    /// Throws if not found or init/ctor not yet called.
    /// @param variable
    /// @return
    T& operator[](const std::string& variable);

    /// Determine if a variable is in the storage. Uses _s for compile time hash
    /// @param hash
    /// @return
    bool  has(const uint64_t& hash);
    /// Determine if a variable is in the storage.
    /// @param variable
    /// @return
    bool has(const std::string& variable);

    /// Returns a list of the variables stored
    /// @return
    std::vector<std::string> variables();

    /// Initialize the storage with a set of variables. Values default to -9999  or nullptr
    /// Not defined for unsigned values
    /// @param variables
    void init(std::set<std::string>& variables);

    /// Returns the number of variables stored
    /// @return
    size_t size();

  private:
    const uint64_t seed = 2654435761U;
    template <typename Item> class wyandFunctor
    {
      public:
        uint64_t operator ()  (const Item& key, uint64_t seed = 2654435761U) const
        {
            return wyhash(&key, sizeof(Item), seed);
        }

    };
    typedef wyandFunctor<uint64_t> hasher_t;
    typedef boomphf::mphf< uint64_t, hasher_t  > boophf_t;

    // Holds the name-value pair in the variable store hashmap
    // we do this as we hold a hash and not the name
    struct var
    {
        T value{};
        double xxhash; // holds the xxhash value so we can confirm we get the right thing back from BBHash
        std::string variable;
    };
    // Note that we have to explicitly check if what we get back is what we wanted as
    // mphf do not guarantee what asking for something outside of the map returns a sane answer
    // https://github.com/rizkg/BBHash/issues/12

    // perfect hashfn + variable storage
    std::unique_ptr<boophf_t> _variable_bphf;
    std::vector<var> _variables;

    // Total number of variables stored
    size_t _size;

};


template<typename T>
variablestorage<T>::variablestorage()
{
    _size = 0;
    _variable_bphf = nullptr;
}

template<typename T>
variablestorage<T>::~variablestorage()
{

}

template<typename T>
T& variablestorage<T>::operator[](const uint64_t& hash)
{
    if(!_variable_bphf)
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("Variable " + std::to_string(hash) + " does not exist."));
    }
    uint64_t  idx = _variable_bphf->lookup(hash);

    // did the table return garabage?
    //mphf might return an index, but it isn't actually what we want. double check the hash
    if (idx >=  _size ||
        _variables[idx].xxhash != hash )
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("Variable " + std::to_string(hash) + " does not exist."));

    return _variables[idx].value;
}

template<typename T>
T& variablestorage<T>::operator[](const std::string& variable)
{
    if(!_variable_bphf)
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("Variable " + variable + " does not exist."));
    }

    uint64_t hash = xxh64::hash (variable.c_str(), variable.length(), this->seed);
    uint64_t  idx = _variable_bphf->lookup(hash);

    // did the table return garabage?
    //mphf might return an index, but it isn't actually what we want. double check the hash
    if( idx >=  _size ||
        _variables[idx].xxhash != hash)
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("Variable " + variable + " does not exist."));

    return _variables[idx].value;
}

template<typename T>
bool variablestorage<T>::has(const uint64_t& hash)
{
    if (_size == 0) return false;

    uint64_t  idx = _variable_bphf->lookup(hash);

    // did the table return garabage?
    //mphf might return an index, but it isn't actually what we want. double check the hash
    if(_variables[idx].xxhash != hash)
        return false;

    return true;
}

template<typename T>
bool variablestorage<T>::has(const std::string& variable)
{
    uint64_t hash = xxh64::hash (variable.c_str(), variable.length(), 2654435761U);
    return has(hash);

};

template<typename T>
std::vector<std::string> variablestorage<T>::variables()
{
    std::vector<std::string> vars;
    for(auto itr:_variables)
    {
        vars.push_back(itr.variable);
    }
    return vars;

}

template<typename T>
variablestorage<T>::variablestorage(std::set<std::string>& variables)
    : variablestorage()
{

    init(variables);
}

template<typename T>
void variablestorage<T>::init(std::set<std::string>& variables)
{

    std::vector<u_int64_t> hash_vec;
    for(auto& v : variables)
    {
        uint64_t hash = xxh64::hash (v.c_str(), v.length(), seed);
        hash_vec.push_back(hash);
    }

    _variable_bphf = std::unique_ptr<boomphf::mphf<u_int64_t,hasher_t>>(
        new boomphf::mphf<u_int64_t,hasher_t>(hash_vec.size(),hash_vec,1,2,false,false));

    _variables.resize(variables.size());
    for(auto& v : variables)
    {
        uint64_t hash = xxh64::hash (v.c_str(), v.length(), seed);
        uint64_t  idx = _variable_bphf->lookup(hash);

        _variables[idx].variable = v;
        _variables[idx].xxhash = hash;
    }

    _size = variables.size();
}


template<typename T>
size_t variablestorage<T>::size()
{
    return _size;
}
