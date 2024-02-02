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

// Needed for the SPDLOG_ macro calls
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

//spdlog includes
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

 enum log_level
 {
     verbose,
     debug,
     warning,
     info,
     error
 };

//#ifdef USE_MPI
//#define MPI_RANK_DBG(RANK) if(_comm_world.rank() == RANK) \
//    {\
//    LOG_DEBUG << "\n\n-----------------------------\n\nI am PID " << getpid() <<", attach within 45s\n\n-----------------------------\n\n"; \
//            sleep(45);\
//    }
//#endif