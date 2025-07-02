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
#include <vector>
#include <string>
#include <stdexcept>
#include <cstring>

class ParallelUGRIDWriter {
public:
    ParallelUGRIDWriter(const std::string& filename,
                        boost::mpi::communicator& comm,
                        int nNodes, int nFaces, int nTimeSteps);

    void write_mesh_topology(const std::vector<double>& node_x,
                             const std::vector<double>& node_y,
                             const std::vector<double>& node_z,
                             const std::vector<int>& face_nodes,
                             const std::vector<int>& face_neighbors,
                             const std::vector<int>& face_global_id,
                             const std::vector<double>& static_param);

    void write_time_step(int time_idx,
                         const std::vector<double>& var_data,
                         const std::string& var_name);

    void close();

    ~ParallelUGRIDWriter();

private:
    int ncid;
    boost::mpi::communicator& comm;
    int rank, size;
    int nNodes, nFaces, nTimeSteps;

    int node_dimid, face_dimid, time_dimid, face_corner_dimid;
    int node_x_varid, node_y_varid, node_z_varid;
    int face_nodes_varid, face_neighbors_varid, face_gid_varid;
    int static_param_varid;
    std::vector<int> time_varids;
    std::vector<std::string> time_varnames;

    void define_dimensions();
    void define_variables();
    void define_attributes();
    void check_nc(int stat, const std::string& msg);
};

