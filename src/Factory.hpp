////////////////////////////////////////////////////////////////////////////
// Factory.hop

#pragma once

// #include <memory>
#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <functional>

template <class Interface, class... ConstructorArgs>
class Factory {
public:
  // Only ever hand out unique pointers
  static boost::shared_ptr<Interface> Create(std::string name, ConstructorArgs... args);
  // Actual registration function
  static void RegisterFactoryFunction(std::string name,
        std::function<Interface*(ConstructorArgs...)> implementationFactoryFunction);
  // Get all keys from the registry
  static std::vector<std::string> registered_keys(void);

private:
  // the actual registry is private to this class
  static std::map<std::string, std::function<Interface*(ConstructorArgs...)>>& registry();

  Factory(){};
  // Do NOT implement copy methods
  Factory(Factory const& copy);
  Factory& operator=(Factory const& copy);
};
// Registrar helper class for Factory
template <class Interface, class Implementation, class... ConstructorArgs>
class Registrar {
public:
  Registrar(std::string className)
  {
    Factory<Interface,ConstructorArgs...>::RegisterFactoryFunction(className, [](ConstructorArgs... args) -> Interface * {return new Implementation(args...);});
  }
};

template <class Interface, class... ConstructorArgs>
boost::shared_ptr<Interface> Factory<Interface,ConstructorArgs...>::Create(std::string name, ConstructorArgs... args)
{
  Interface * instance = nullptr;

  // find the name in the registry and call factory method.
  auto it = registry().find(name);
  if(it != registry().end())
    instance = it->second(args...);

  // wrap instance in a unique ptr and return (if created)
  if(instance == nullptr)
    throw "Table name not found in registry."; // TODO better exception
  return boost::shared_ptr<Interface>(instance);
}

template <class Interface, class... ConstructorArgs>
std::vector<std::string> Factory<Interface,ConstructorArgs...>::registered_keys()
{
  // copy all keys from the registry map into a vector
  std::vector<std::string> keys;
  for (auto const& elem : registry() ) {
    keys.push_back(elem.first);
  }
  return keys;
}

template <class Interface, class... ConstructorArgs>
void Factory<Interface,ConstructorArgs...>::RegisterFactoryFunction(std::string name,
std::function<Interface*(ConstructorArgs...)> implementationFactoryFunction)
{
  // register a derived class factory function
  registry()[name] = implementationFactoryFunction;
}

template <class Interface, class... ConstructorArgs>
std::map<std::string, std::function<Interface*(ConstructorArgs...)> >& Factory<Interface,ConstructorArgs...>::registry()
{
  // Get the singleton instance of the registry map
  static std::map<std::string, std::function<Interface*(ConstructorArgs...)> > registry;
  return registry;
}
