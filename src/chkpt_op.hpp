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

#include <boost/filesystem/path.hpp>
#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>

#include "hpc_scheduler_info.hpp"
#include "logger.hpp"
#include "timeseries/netcdf.hpp"


// Checkpointing options
class chkptOp
{
  public:
    chkptOp():
                do_checkpoint{false},
                load_from_checkpoint{false},
                on_last{false},
                checkpoint_request_terminate{false}
    {
        abort_when_wallclock_left = boost::posix_time::minutes(2);
    }

    boost::filesystem::path ckpt_path; // root path to chckpoint folder
    netcdf in_savestate; // if we are loading from checkpoint
    bool do_checkpoint; // should we check point?
    bool load_from_checkpoint; // are we loading from a checkpoint?
    // amount of time to give ourselves to bail and checkpoint if we have a wall clock limit
    boost::posix_time::time_duration abort_when_wallclock_left;
    boost::optional<bool> on_outta_time; // bail when we are out of time
    boost::optional<bool> on_last; //only checkpoint on the last timestep
    boost::optional<size_t> frequency; // frequency of checkpoints


    boost::optional<boost::posix_time::ptime> specific_datetime; // at a specific date-time
    boost::optional<boost::posix_time::ptime> specific_time; // at a specific time

    // used to stop the simulation when we checkpoint when we are outta time
    bool checkpoint_request_terminate;

    /**
     * Should checkpointing occur
     * @param current_ts
     * @param is_last_ts
     * @return
     */
    bool should_checkpoint(size_t current_ts,
        bool is_last_ts,
        hpc_scheduler_info& scheduler_info,
        boost::mpi::communicator& comm_world,
        const boost::posix_time::ptime& _current_date)
    {
        if(!do_checkpoint)
            return false;

        if(on_last && *on_last && is_last_ts)
            return true;

        // don't checkpoint on the first ts if we are doing frequency checkpoints
        if( frequency && current_ts !=0 && (current_ts % *frequency ==0) )
            return true;

        // check if we are running out of time
        if(on_outta_time && *on_outta_time &&
            scheduler_info.has_wallclock_limit 
            )
        {
            
            int outoftime =  scheduler_info.wallclock_remaining() <= abort_when_wallclock_left;
            int global_outoftime = -1;

            // find is anyone thinks we should bail
            boost::mpi::all_reduce(comm_world, outoftime, global_outoftime, boost::mpi::maximum<int>());

            if(global_outoftime)
            {
                SPDLOG_DEBUG("Detected wallclock of {} remaining. Triggering checkpoint.",
                            boost::posix_time::to_simple_string(scheduler_info.wallclock_remaining()));
                checkpoint_request_terminate = true;
                return true;
            }

        }

        if(specific_time)
        {
            if( (_current_date.time_of_day().hours() == specific_time->time_of_day().hours()) &&
                (_current_date.time_of_day().minutes() == specific_time->time_of_day().minutes()) )
                return true;
        }

        if(specific_datetime)
        {
            if(_current_date == *specific_datetime)
                return true;
        }

        return false;
    }
};
