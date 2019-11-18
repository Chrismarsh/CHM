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

//std includes
#include <string>
#include <set>
#include <vector>

//boost includes
#include <boost/date_time/posix_time/posix_time.hpp> // for boost::posix

//Gdal includes
#include <ogr_spatialref.h>

//CHM includes
#include "exception.hpp"
#include "logger.hpp"
#include "station.hpp"
#include "netcdf.hpp"
#include "timeseries.hpp"
#include "triangulation.hpp"
#include "filter_base.hpp"
/**
 * Main meteorological data coordinator. Opens from a variety of sources and ensures that each virtual station has this timestep's information
 * regardless of the source data type.
 */
class metdata
{
  public:

    /// Holds the metadata for an ascii file to load
    ///
    struct ascii_metdata
    {
        double latitude, longitude, elevation;
        std::string path;
        std::string id;


        //if we use text file inputs, each station can have its own filer (ie., winds at different heights). So we need to save the filter
        //and run it on a per-station config.
        std::vector<boost::shared_ptr<filter_base>> filters;
    };

    struct ascii_data
    {
        //if we use text file inputs, each station can have its own filer (ie., winds at different heights). So we need to save the filter
        //and run it on a per-station config.
        std::vector<boost::shared_ptr<filter_base>> filters;

        std::string id;
        // these are loaded into by metdata. Essentially this becomes like the old station
        timeseries _obs;
        timeseries::iterator _itr;
    };

    metdata();
    metdata(boost::shared_ptr< triangulation > mesh);
    ~metdata();

    /// Loads a netcdf file. Must be a 2D structured grid of stations. Expects times to be in UTC+0
    /// @param path
    /// @param filters
    void load_from_netcdf(const std::string& path, std::map<std::string, boost::shared_ptr<filter_base> > filters);

    /// Loads the standard ascii timeseries. Needs to be in UTC+0
    /// @param path
    /// @param filters
    /// @param utc_offset Positive offset going west. So the normal UTC-6 would be UTC_offset:6
    void load_from_ascii(std::vector<ascii_metdata> stations, int utc_offset);

    void write_stations_to_ptv(const std::string& path);

    /// For all the stations loaded from ascii files, find the latest start time, and the earliest end time that is consistent across all stations
    /// @return
    std::pair<boost::posix_time::ptime, boost::posix_time::ptime> start_end_times();

    /// Number of stations
    /// @return
    size_t nstations();

    std::shared_ptr<station> at(size_t idx);

    boost::posix_time::ptime start_time();
    boost::posix_time::ptime end_time();

    /// Subsets all timeseries to begin at [start, end]. For ascii, the underlying timeseries is modified.
    /// For nc, internal offsets are computed to start, end.
    /// This updates the internal start and end times, as well as resets the current time to be = start
    /// @param start
    /// @param end
    void subset(boost::posix_time::ptime start, boost::posix_time::ptime end);

    /// Check that all ascii stations have the same start/end times
    /// @return
    void check_ts_consistency();

    /// Timestep duration in seconds
    /// @return
    size_t dt();

    /// Populates the stations' with the next timesteps' value
    /// @return False if no more timesteps
    bool next();

    /// Advances 1 timestep in the netcdf files
    bool next_nc();

    /// Advances 1 timestep from the ascii timeseries
    bool next_ascii();

  private:

    // NetCDF specific variables
    // -----------------------------------
        //if we use netcdf, store it here
        std::unique_ptr<netcdf> _nc;

        //if we use netcdf, we need to save the filters and run it once every timestep.
        std::map<std::string, boost::shared_ptr<filter_base>>_netcdf_filters;

        std::set<std::string> _provides_from_nc_filters;

        // if false, we are using ascii files
        bool _use_netcdf;

    // -----------------------------------
    // ASCII met data specific variables

        //Essentially what the old stations turned into.
        // Holds all the met data to init that stations + the underlying timeseries data
        // Mapped w/ stations ID -> metdata
        std::map<std::string, std::unique_ptr<ascii_data>> _ascii_stations;




    // -----------------------------------

    // This is a different approach than how stations used to work
    // Now, they only hold the current timestep, which is refilled every model timestep by metdata
    // Our sources of data are netcdf, or ascii (or whatever in the future)
    std::vector< std::shared_ptr<station>> _stations;


    // Total number of stations
    size_t _nstations;

    //ptr to the core:: owned mesh. Metdata needs it to know what it should load when in MPI mode
    boost::shared_ptr< triangulation > _mesh;

    boost::posix_time::ptime _start_time, _end_time;
    boost::posix_time::ptime _current_ts;
    boost::posix_time::time_duration _dt;

    // computes the dt
    void compute_dt();

};

