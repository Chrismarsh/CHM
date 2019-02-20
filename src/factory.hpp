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

#include <map>
#include <vector>
#include <functional>

#include <boost/shared_ptr.hpp>

#include "exception.hpp"

template <class Interface, class... ConstructorArgs>
class factory {
public:
  // Only ever hand out unique pointers
  static boost::shared_ptr<Interface> create(std::string name, ConstructorArgs... args);
  // Actual registration function
  static void register_factory_function(std::string name,
        std::function<Interface*(ConstructorArgs...)> implementation_constructor_function);
  // Get all keys from the registry
  static std::vector<std::string> registered_keys(void);

private:
  // the actual registry is private to this class
  static std::map<std::string, std::function<Interface*(ConstructorArgs...)>>& registry();

  factory(){};
  // Remove copy constructor methods
  factory(factory const& copy) = delete;
  factory& operator=(factory const& copy) = delete;
};
// registrar helper class for factory
template <class Interface, class Implementation, class... ConstructorArgs>
class registration_helper {
public:
  registration_helper(std::string className)
  {
    factory<Interface,ConstructorArgs...>::register_factory_function(className, [](ConstructorArgs... args) -> Interface * {return new Implementation(args...);});
  }
};

template <class Interface, class... ConstructorArgs>
boost::shared_ptr<Interface> factory<Interface,ConstructorArgs...>::create(std::string name, ConstructorArgs... args)
{
  Interface * instance = nullptr;

  // find the name in the registry and call factory method.
  auto it = registry().find(name);
  if(it != registry().end())
    instance = it->second(args...);

  // wrap instance in a shared ptr and return (if created)
  if(instance == nullptr) {
    CHM_THROW_EXCEPTION(module_not_found, "Requested module not found in registry: [" + name +"]");
  }
  return boost::shared_ptr<Interface>(instance);
}

template <class Interface, class... ConstructorArgs>
std::vector<std::string> factory<Interface,ConstructorArgs...>::registered_keys()
{
  // copy all keys from the registry map into a vector
  std::vector<std::string> keys;
  for (auto const& elem : registry() ) {
    keys.push_back(elem.first);
  }
  return keys;
}

template <class Interface, class... ConstructorArgs>
void factory<Interface,ConstructorArgs...>::register_factory_function(std::string name,
std::function<Interface*(ConstructorArgs...)> implementation_constructor_function)
{
  // register a derived class factory function
  registry()[name] = implementation_constructor_function;
}

template <class Interface, class... ConstructorArgs>
std::map<std::string, std::function<Interface*(ConstructorArgs...)> >& factory<Interface,ConstructorArgs...>::registry()
{
  // Get the singleton instance of the registry map
  static std::map<std::string, std::function<Interface*(ConstructorArgs...)> > registry;
  return registry;
}
