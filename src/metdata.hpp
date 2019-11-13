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
    struct ascii_data
    {
        double latitude, longitude, elevation;
        std::string path;
        std::string name;
        std::vector<boost::shared_ptr<filter_base>> filters;
    };

    metdata();
    metdata(const boost::shared_ptr< triangulation > mesh);
    ~metdata();

    /// Loads a netcdf file. Must be a 2D structured grid of stations. Expects times to be in UTC+0
    /// @param path
    /// @param filters
    void load_from_netcdf(const std::string& path, std::map<std::string, boost::shared_ptr<filter_base> > filters);

    /// Loads the standard ascii timeseries. Needs to be in UTC+0
    /// @param path
    /// @param filters
    /// @param utc_offset Positive offset going west. So the normal UTC-6 would be UTC_offset:6
    void load_from_ascii(std::vector<ascii_data> stations, int utc_offset);

    void write_stations_to_ptv(const std::string& path);

    /// Number of stations
    /// @return
    size_t nstations();

    std::shared_ptr<station> at(size_t idx);
  private:

    // NetCDF specific variables

    //if we use netcdf, store it here
    std::unique_ptr<netcdf> _nc;

    //if we use netcdf, we need to save the filters and run it once every timestep.
    std::map<std::string, boost::shared_ptr<filter_base>>_netcdf_filters;

    std::set<std::string> _provides_from_nc_filters;

    // if false, we are using ascii files
    bool _use_netcdf;

    //if we use ascii files, store them here
    std::vector< std::unique_ptr<timeseries> > _ascii_timeseries;

    std::vector< std::shared_ptr<station>> _stations;



    // Total number of stations
    size_t _nstations;

    //ptr to the core:: owned mesh. Metdata needs it to know what it should load when in MPI mode
    const boost::shared_ptr< triangulation > _mesh;



    //if we use text file inputs, each station can have its own filer (ie., winds at different heights). So we need to save the filter
    //and run it on a per-station config.
    std::map<std::string, std::vector<boost::shared_ptr<filter_base>> > _txtmet_filters;

};

