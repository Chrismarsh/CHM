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

#include "variablestorage.hpp"

variablestorage::variablestorage()
{
    _size = 0;
    _variable_bphf = nullptr;
}
variablestorage::~variablestorage()
{

}

double& variablestorage::operator[](const uint64_t& hash)
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


double& variablestorage::operator[](const std::string& variable)
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

bool variablestorage::has(const uint64_t& hash)
{
    if (_size == 0) return false;

    uint64_t  idx = _variable_bphf->lookup(hash);

    // did the table return garabage?
    //mphf might return an index, but it isn't actually what we want. double check the hash
    if(_variables[idx].xxhash != hash)
        return false;

    return true;
}

bool variablestorage::has(const std::string& variable)
{
    uint64_t hash = xxh64::hash (variable.c_str(), variable.length(), 2654435761U);
    return has(hash);

};

std::vector<std::string> variablestorage::variables()
{
    std::vector<std::string> vars;
    for(auto itr:_variables)
    {
        vars.push_back(itr.variable);
    }
    return vars;

}

variablestorage::variablestorage(std::set<std::string>& variables)
: variablestorage()
{

    init(variables);
}
void variablestorage::init(std::set<std::string>& variables)
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
    for(auto v : variables)
    {
        uint64_t hash = xxh64::hash (v.c_str(), v.length(), seed);
        uint64_t  idx = _variable_bphf->lookup(hash);
        _variables[idx].value = -9999.0;
        _variables[idx].variable = v;
        _variables[idx].xxhash = hash;
    }

    _size = variables.size();
}
size_t variablestorage::size()
{
    return _size;
}