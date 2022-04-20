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



#include <string>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/ptree.hpp>

namespace pt = boost::property_tree;

#include <tbb/concurrent_vector.h>

#include "interpolation.hpp"

#include "math/coordinates.hpp"


/**
 * Basin wide parameters such as transmissivity, solar elevation, solar aspect, etc.
 *
 **/
class global
{

    struct Distance
    {

    };



    //want to let core modify date time, etc without showing a public interface.
    //This is because global gets passed to all modules and a rogue module could do something dumb
    //const doesn't save us as we actually do want to modify things
    friend class core;

private:
    boost::posix_time::ptime _current_date;

    int _dt; //seconds
    bool _is_geographic;
    bool _is_point_mode;
    bool _from_checkpoint;


public:

    bool is_point_mode();

    // UTC offset
    int _utc_offset;
    bool is_geographic();
    global();
    int year();
    int day();
    interp_alg interp_algorithm;
    /*
     * Month on [1,12]
     */
    int month();
    int hour();
    int min();
    int sec();
    int dt();
    boost::posix_time::ptime posix_time();
    uint64_t posix_time_int();

    size_t timestep_counter; // the timestep we are on, start = 0


    bool first_time_step;

    /**
     * Returns true if we have loaded from a check point file
     * @return
     */
    bool from_checkpoint();


    pt::ptree parameters;


};
