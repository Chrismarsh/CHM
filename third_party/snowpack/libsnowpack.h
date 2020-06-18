/*
 *  SNOWPACK stand-alone
 *
 *  Copyright WSL Institute for Snow and Avalanche Research SLF, DAVOS, SWITZERLAND
*/
/*  This file is part of libsnowpack.
    libsnowpack is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libsnowpack is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with libsnowpack.  If not, see <http://www.gnu.org/licenses/>.
*/
/**
 * @file libsnowpack.h
 * @version 10.02
 * This is the header file to include for projects making use of libsnowpack
 */

#ifndef LIBSNOWPACK_H
#define LIBSNOWPACK_H

#include "Constants.h"
#include "DataClasses.h"
#include "Hazard.h"
#include "Laws_sn.h"
#include "Meteo.h"
#include "Saltation.h"
#include "SnowDrift.h"
#include "SnowpackConfig.h"
#include "Stability.h"
#include "Utils.h"

#include "plugins/AsciiIO.h" //for direct calls to AsciiIO
#include "plugins/SmetIO.h"  //for direct calls to SmetIO
#include "plugins/SnowpackIO.h"
#include "plugins/SnowpackIOInterface.h"

#include "snowpackCore/Aggregate.h"
#include "snowpackCore/Canopy.h"
#include "snowpackCore/Metamorphism.h"
#include "snowpackCore/PhaseChange.h"
#include "snowpackCore/ReSolver1d.h"
#include "snowpackCore/Snowpack.h"
#include "snowpackCore/Solver.h"
#include "snowpackCore/WaterTransport.h"

#endif

