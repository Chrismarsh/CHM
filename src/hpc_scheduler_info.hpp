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

#include <cstdlib>
#include <string>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/posix_time_io.hpp>

#include "exception.hpp"
#include "logger.hpp"


// Detects various information about the HPC scheduler we might be running under.
class hpc_scheduler_info
{
public:

  boost::posix_time::time_duration max_wallclock; // maximum wallclock in seconds
  boost::posix_time::ptime wallclock_start; // time we started the simulation at
  bool has_wallclock_limit;
  std::string job_name;

    hpc_scheduler_info()
    {
        max_wallclock = boost::posix_time::seconds(0);
        has_wallclock_limit = false;
        job_name = "";
    }

    /**
     * If we have a wallclock limit, how much time left?
     * Only produces a useful delta if has_wallclock_limit = true;
     * @return
     */
    boost::posix_time::time_duration wallclock_remaining()
    {
        return max_wallclock - (boost::posix_time::second_clock::local_time()-wallclock_start);
    }

    void detect()
    {
        // Check if we are running under slurm
        const char* SLURM_JOB_ID = std::getenv("SLURM_JOB_ID");
        if (SLURM_JOB_ID)
        {
            const char* SLURM_TASK_PID = std::getenv("SLURM_TASK_PID"); // The process ID of the task being started.
            const char* SLURM_PROCID =
                std::getenv("SLURM_PROCID"); // The MPI rank (or relative process ID) of the current process
            job_name = SLURM_JOB_ID;
            SPDLOG_DEBUG("Detected running under SLURM as jobid {}", job_name);
            SPDLOG_DEBUG("SLURM_TASK_PID = {}", SLURM_TASK_PID);
            SPDLOG_DEBUG("SLURM_PROCID = {} ", SLURM_PROCID);
        }


        // check if we are running under PBS
        const char* PBS_JOB_ID = std::getenv("PBS_JOBID");
        if(PBS_JOB_ID)
        {
            job_name = PBS_JOB_ID;
            SPDLOG_DEBUG("Detected running under PBS as jobid {}", job_name);

        }

        const char* CHM_WALLCLOCK = std::getenv("CHM_WALLCLOCK_LIMIT");
        if(CHM_WALLCLOCK)
        {
            try {
                max_wallclock = boost::posix_time::duration_from_string(CHM_WALLCLOCK);
                has_wallclock_limit = true;
                wallclock_start = boost::posix_time::second_clock::local_time();
                SPDLOG_DEBUG("Detected a max wallclock of {}", boost::posix_time::to_simple_string(max_wallclock));
            } catch (...) {
                CHM_THROW_EXCEPTION(chm_error, "The value given for environment variable CHM_WALLCLOCK is invalid");
            }
        }
    }
};
