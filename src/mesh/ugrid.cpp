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


#include "ugrid.hpp"

ParallelUGRIDWriter::ParallelUGRIDWriter(const std::string& filename,
                                         boost::mpi::communicator& comm_,
                                         int nNodes_, int nFaces_, int nTimeSteps_)
    : comm(comm_), nNodes(nNodes_), nFaces(nFaces_), nTimeSteps(nTimeSteps_) {
    rank = comm.rank();
    size = comm.size();

    // Open NetCDF file in parallel
    check_nc(nc_create_par(filename.c_str(), NC_NETCDF4 | NC_CLOBBER, comm_, MPI_INFO_NULL, &ncid),
             "nc_create_par");

    define_dimensions();
    define_variables();
    define_attributes();

    check_nc(nc_enddef(ncid), "nc_enddef");
}

void ParallelUGRIDWriter::define_dimensions() {
    check_nc(nc_def_dim(ncid, "nNodes", nNodes, &node_dimid), "def_dim nNodes");
    check_nc(nc_def_dim(ncid, "nFaces", nFaces, &face_dimid), "def_dim nFaces");
    check_nc(nc_def_dim(ncid, "nTime", nTimeSteps, &time_dimid), "def_dim nTime");
    check_nc(nc_def_dim(ncid, "nFaceCorners", 3, &face_corner_dimid), "def_dim nFaceCorners");
}

void ParallelUGRIDWriter::define_variables() {
    int dim1[1] = {node_dimid};
    check_nc(nc_def_var(ncid, "node_x", NC_DOUBLE, 1, dim1, &node_x_varid), "def_var node_x");
    check_nc(nc_def_var(ncid, "node_y", NC_DOUBLE, 1, dim1, &node_y_varid), "def_var node_y");
    check_nc(nc_def_var(ncid, "node_z", NC_DOUBLE, 1, dim1, &node_z_varid), "def_var node_z");

    int dim2[2] = {face_dimid, face_corner_dimid};
    check_nc(nc_def_var(ncid, "face_nodes", NC_INT, 2, dim2, &face_nodes_varid), "def_var face_nodes");
    check_nc(nc_def_var(ncid, "face_neighbors", NC_INT, 2, dim2, &face_neighbors_varid), "def_var face_neighbors");

    int dimf[1] = {face_dimid};
    check_nc(nc_def_var(ncid, "face_global_id", NC_INT, 1, dimf, &face_gid_varid), "def_var face_global_id");
    check_nc(nc_def_var(ncid, "static_param", NC_DOUBLE, 1, dimf, &static_param_varid), "def_var static_param");
}

void ParallelUGRIDWriter::define_attributes() {
    int mesh_varid;
    check_nc(nc_def_var(ncid, "Mesh2", NC_INT, 0, nullptr, &mesh_varid), "def_var Mesh2");
    check_nc(nc_put_att_text(ncid, mesh_varid, "cf_role", strlen("mesh_topology"), "mesh_topology"), "att cf_role");
    int topo_dim = 2;
    check_nc(nc_put_att_int(ncid, mesh_varid, "topology_dimension", NC_INT, 1, &topo_dim), "att topology_dimension");
    check_nc(nc_put_att_text(ncid, mesh_varid, "node_coordinates", strlen("node_x node_y node_z"), "node_x node_y node_z"), "att node_coordinates");
    check_nc(nc_put_att_text(ncid, mesh_varid, "face_node_connectivity", strlen("face_nodes"), "face_nodes"), "att face_node_connectivity");
    check_nc(nc_put_att_text(ncid, mesh_varid, "face_dimension", strlen("nFaces"), "nFaces"), "att face_dimension");
    int start_index = 0;
    check_nc(nc_put_att_int(ncid, mesh_varid, "start_index", NC_INT, 1, &start_index), "att start_index");
}

void ParallelUGRIDWriter::write_mesh_topology(const std::vector<double>& node_x,
                                              const std::vector<double>& node_y,
                                              const std::vector<double>& node_z,
                                              const std::vector<int>& face_nodes,
                                              const std::vector<int>& face_neighbors,
                                              const std::vector<int>& face_global_id,
                                              const std::vector<double>& static_param) {
    // All ranks must enter data mode
    check_nc(nc_var_par_access(ncid, node_x_varid, NC_COLLECTIVE), "par_access node_x");
    check_nc(nc_var_par_access(ncid, node_y_varid, NC_COLLECTIVE), "par_access node_y");
    check_nc(nc_var_par_access(ncid, node_z_varid, NC_COLLECTIVE), "par_access node_z");
    check_nc(nc_var_par_access(ncid, face_nodes_varid, NC_COLLECTIVE), "par_access face_nodes");
    check_nc(nc_var_par_access(ncid, face_neighbors_varid, NC_COLLECTIVE), "par_access face_neighbors");
    check_nc(nc_var_par_access(ncid, face_gid_varid, NC_COLLECTIVE), "par_access face_gid");
    check_nc(nc_var_par_access(ncid, static_param_varid, NC_COLLECTIVE), "par_access static_param");

    // Write node coordinates (all nodes must be present on all ranks)
    check_nc(nc_put_var_double(ncid, node_x_varid, node_x.data()), "put_var node_x");
    check_nc(nc_put_var_double(ncid, node_y_varid, node_y.data()), "put_var node_y");
    check_nc(nc_put_var_double(ncid, node_z_varid, node_z.data()), "put_var node_z");

    // Write face data (each rank writes its local faces)
    size_t start_face = 0; // For a real distributed mesh, set this to the global offset for this rank
    size_t count_face = face_global_id.size();

    check_nc(nc_put_vara_int(ncid, face_nodes_varid, &start_face, &count_face, face_nodes.data()), "put_vara face_nodes");
    check_nc(nc_put_vara_int(ncid, face_neighbors_varid, &start_face, &count_face, face_neighbors.data()), "put_vara face_neighbors");
    check_nc(nc_put_vara_int(ncid, face_gid_varid, &start_face, &count_face, face_global_id.data()), "put_vara face_gid");
    check_nc(nc_put_vara_double(ncid, static_param_varid, &start_face, &count_face, static_param.data()), "put_vara static_param");
}

void ParallelUGRIDWriter::write_time_step(int time_idx,
                                          const std::vector<double>& var_data,
                                          const std::string& var_name) {
    int varid;
    int dimids[2] = {time_dimid, face_dimid};
    // If variable not yet defined, define it
    if (std::find(time_varnames.begin(), time_varnames.end(), var_name) == time_varnames.end()) {
        check_nc(nc_redef(ncid), "redef for time var");
        check_nc(nc_def_var(ncid, var_name.c_str(), NC_DOUBLE, 2, dimids, &varid), "def_var time_var");
        check_nc(nc_enddef(ncid), "enddef for time var");
        time_varids.push_back(varid);
        time_varnames.push_back(var_name);
    } else {
        varid = time_varids[std::distance(time_varnames.begin(),
                                          std::find(time_varnames.begin(), time_varnames.end(), var_name))];
    }
    // Write data for this time step (each rank writes its local faces)
    size_t start[2] = {static_cast<size_t>(time_idx), 0}; // 0: replace with global face offset for this rank
    size_t count[2] = {1, var_data.size()};
    check_nc(nc_put_vara_double(ncid, varid, start, count, var_data.data()), "put_vara time_var");
}

void ParallelUGRIDWriter::close() {
    if (ncid >= 0)
        nc_close(ncid);
}

ParallelUGRIDWriter::~ParallelUGRIDWriter() {
    close();
}

void ParallelUGRIDWriter::check_nc(int stat, const std::string& msg) {
    if (stat != NC_NOERR) {
        throw std::runtime_error(msg + ": " + nc_strerror(stat));
    }
}
