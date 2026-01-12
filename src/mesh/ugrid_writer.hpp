// 	Copyright (C) 2011  Chris Marsh
//
// 	This program is free software: you can redistribute it and/or modify
// 	it under the terms of the GNU General Public License as published by
// 	the Free Software Foundation, either version 3 of the License, or
// 	(at your option) any later version.
//
// 	This program is distributed in the hope that it will be useful,
// 	but WITHOUT ANY WARRANTY; without even the implied warranty of
// 	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// 	GNU General Public License for more details.
//
// 	You should have received a copy of the GNU General Public License
// 	along with this program.  If not, see <http://www.gnu.org/licenses/>.

#pragma once

#include <netcdf.h>
#include <netcdf_par.h>
#include <boost/mpi.hpp>
#include <memory>
#include <algorithm>
#include <vector>
#include <string>
#include <utility>
#include <cstdint>
#include <format>

#include <boost/mpi.hpp>
#include <boost/optional.hpp>

#include "triangulation.hpp"
#include "timer.hpp"
#include "global.hpp"


class ugrid_writer
{
public:
    ugrid_writer(mesh m, std::shared_ptr<global> g, bool write_parameters, std::string fname, bool use_zarr = false);
    ~ugrid_writer();

    /**
     * Check a netcdf C call's return value and convert to exception if needed
     * @param status
     */
    void nc_chk_ret(int status);

    void close_ugrid();
    void write_ugrid(const std::vector<std::string>& output_variables);
    void open_ugrid(const std::vector<std::string>& output_variables);
    void init_ugrid(const std::vector<std::string>& output_variables);
    void set_output_cadence(const boost::optional<size_t>& frequency,
                            const boost::optional<size_t>& only_last_n,
                            const boost::optional<size_t>& rotate_frequency);
    void set_chunking_override(const boost::optional<size_t>& chunk_len_steps,
                               const boost::optional<double>& chunk_target_mb);
    void set_store_path(std::string store_path);
    const std::string& store_path() const;
    // Probe the existing ugrid file to determine how many timesteps are already written.
    size_t probe_time_index();

    bool bitgroom;
    bool compress;
private:
    std::string build_store_uri(const std::string& store_path) const;
    bool store_exists() const;
    size_t compute_time_chunk_len(size_t max_faces_per_rank, size_t num_output_vars) const;
    size_t read_time_index(int fid) const;

    //holds the file id for the ugrid output netcdf
    int _ugrid_fid;

    //maps the variable string to the netcdf id to write to file
    std::map<std::string, int> _ugrid_id_var;

    std::string _fname;
    std::string _store_path;
    bool _use_zarr;

    // track the number of outputs that have been done to correctly compute the offset in the nc
    size_t _time_index;

    mesh _mesh;
    std::shared_ptr<global> _global;

    boost::mpi::communicator _comm_world;


    bool _write_parameters;
    boost::optional<size_t> _frequency;
    boost::optional<size_t> _only_last_n;
    boost::optional<size_t> _rotate_frequency;
    boost::optional<size_t> _chunk_len_steps;
    boost::optional<double> _chunk_target_mb;
};
