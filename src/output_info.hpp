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

#include <cstddef>
#include <set>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>
#include <boost/optional.hpp>
#include <memory>
#include <boost/variant.hpp>

#include "logger.hpp"
#include "mesh/ugrid_writer.hpp"
#include "timeseries/timeseries.hpp"
#include "triangulation.hpp"

class vtk_writer;


class output_info
{
public:
    output_info():
    name{""},
    fname{""},
    latitude{0}, longitude{0},
    x{0}, y{0}, type{output_type::output_type_none},
    mesh_output_formats{mesh_outputs::mesh_outputs_none},
    only_last_n{SIZE_MAX},
    write_ghost_neighbors{false}
    {
        face = nullptr;
    }

    enum output_type
    {
        output_type_none,
        time_series,
        mesh
    };
    enum mesh_outputs
    {
        mesh_outputs_none,
        vtp,
        vtu,
        ugrid,
        ascii
    };
    // Should we rorate the ugrid to a new file?
    bool should_rotate(const size_t& max_ts,
                       const size_t& current_ts,
                       const boost::posix_time::ptime& _current_date
                       )
    {
        bool should = false;
        if(rotate_frequency)
        {
            size_t offset = rotate_offset.value_or(0);
            // Use the stored offset from checkpoint so rotation cadence stays aligned even after resume.
            if((current_ts + offset) % *rotate_frequency == 0)
                should = true;
        }

        return should;
    }

    // Should we output?
    bool should_output(const size_t& max_ts,
                       const size_t& current_ts,
                       const boost::posix_time::ptime& _current_date
                       )
    {
        bool should_output = false;

        if(only_last_n)
        {
            auto ts_left = max_ts - current_ts;
            if( ts_left <= *only_last_n) // if we are within the last n timesteps, output
                should_output = true;
        }

        if(frequency)
        {
            if(current_ts % *frequency == 0)
                should_output = true;
        }

        if(specific_time)
        {
            if( (_current_date.time_of_day().hours() == specific_time->time_of_day().hours()) &&
                (_current_date.time_of_day().minutes() == specific_time->time_of_day().minutes()) )
                should_output = true;
        }

        if(specific_datetime)
        {
            if(_current_date == *specific_datetime)
                should_output = true;
        }

        return should_output;

    }

    // print to stdout DEBUG all the valid outputs selected
    void list_outputs()
    {
            SPDLOG_DEBUG("Output frequency options for {}", name);

            if(only_last_n)
                SPDLOG_DEBUG("\tonly_last_n = {}", *only_last_n);
            if(frequency)
                SPDLOG_DEBUG("\tfrequency = {}", *frequency);
            if(specific_time)
                SPDLOG_DEBUG("\tspecific_time = {}", std::to_string(specific_time->time_of_day().hours()) + ":" + std::to_string(specific_time->time_of_day().minutes()));
            if(specific_datetime)
                SPDLOG_DEBUG("\tspecific_datetime = {}", boost::posix_time::to_simple_string(*specific_datetime));
    }

    output_type type; // the type of output, timeseries or mesh
    std::string name;
    mesh_outputs mesh_output_formats;
    std::string fname;
    std::string base_name; // base file name for vtu or ugrid outputs

    // these are input by the user, assumed to be WGS84
    double latitude;
    double longitude;

    // if we are outputting on a projected mesh then we need to store the projected coords here
    double x;
    double y;

    std::set<std::string> variables;
    std::set<std::string> output_parameters;
    mesh_elem face;
    timeseries ts;

    // Output options

    // because the ugrid file can get huge, start a new file every X timesteps
    // defaults to never
    boost::optional<size_t> rotate_frequency;
    // When resuming, keep the original rotation cadence by offsetting current_ts.
    boost::optional<size_t> rotate_offset;

    // every n timesteps
    boost::optional<size_t> frequency;

    // at a specific date-time
    boost::optional<boost::posix_time::ptime> specific_datetime;

    // at a specific time
    boost::optional<boost::posix_time::ptime> specific_time;

    //Only output the last n timesteps. -1 = all
    boost::optional<size_t> only_last_n;

    // include ghost neighbor faces in vtu output
    bool write_ghost_neighbors;

    // bespoke writer to write this output
    // the ugrid contains non-copyable MPI objects so needs to be ptr
    boost::variant< std::shared_ptr<ugrid_writer>, std::shared_ptr<vtk_writer>> writer;

};
