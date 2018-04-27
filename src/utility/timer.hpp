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

#include <chrono>

#include <stack>

//Title: clock
    typedef std::chrono::milliseconds ms;
    typedef std::chrono::nanoseconds ns;
    typedef std::chrono::seconds s;
class timer
{
public:


    /*
    Function: clock
            Creates a clock.

    Parameters: 
    None

    Throws:
    Never

    Returns:   
		
     */
    timer();

    /*
    Function: Tic
            Begins a timer

    Parameters: 
            None

    Throws:
            Never

    Returns:   
				
     */
    void tic();

    /*
    Function: Toc
            Returns the time elapsed in micro-seconds since the last Tic() call

    Parameters: 
            None

    Throws:
            Never

    Returns:   
            double - Time elapsed, micro-seconds. -1.0 if Tic() hasn't been called yet.
     */
template <typename T>
    double toc();

private:
    std::stack<std::chrono::high_resolution_clock::time_point> _timers;
    double _resolution; //in ns
};


template <typename T>
double timer::toc()
{
    auto end = std::chrono::high_resolution_clock::now();
    auto start = _timers.top();
    _timers.pop();
    auto dur = std::chrono::duration_cast< T>(end - start);
    return dur.count();
}