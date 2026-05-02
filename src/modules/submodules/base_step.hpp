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
#include <concepts>

/**
 * @brief Base class for submodules using the Curiously Recurring Template Pattern (CRTP).
 *
 * This class enforces that derived classes implement an `execute_impl(Data&) const` method,
 * and provides a public `execute(Data&) const` method that forwards to the derived implementation.
 *
 * The CRTP pattern allows compile-time polymorphism, enabling static dispatch and
 * avoiding the overhead of virtual functions. This is useful for performance-critical
 * code or when templates are preferred over inheritance.
 *
 * Usage example:
 * @code
 * struct MyStep : public base_step<MyStep, MyData> {
 *     void execute_impl(MyData& d) {
 *         // Step-specific logic here
 *     }
 * };
 * MyStep step;
 * MyData data;
 * step.execute(data); // Calls MyStep::execute_impl(data)
 * @endcode
 */

template<typename T, typename Data>
concept GuaranteeImplementExecute = requires(const T& t, Data& d) {
    { t.execute_impl(d) } -> std::same_as<void>;
};

template<class Derived, class Data>
class base_step {
protected:
    base_step() {};
public:
    void execute(Data& d) const {

        static_assert(GuaranteeImplementExecute<Derived,Data>,
                "Derived class must implement: void execute_impl(Data&) const");

        static_cast<const Derived*>(this)->execute_impl(d);
    }
    ~base_step() {};
};
