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


#include "ugrid_writer.hpp"

ugrid_writer::ugrid_writer(mesh m, boost::shared_ptr<global> g, bool write_parameters, std::string fname):
    _mesh(m),  _fname(fname), _global(g), _write_parameters(write_parameters)
{
    _ugrid_fid = -1;
    _time_index = 0;
    compress = true;
    bitgroom = false;

}

ugrid_writer::~ugrid_writer()
{
    try
    {
        // in theory can throw so log the closure error and swallow the exception
        close_ugrid();
    }catch (const chm_error& e)
    {
        SPDLOG_ERROR(e.what());
    }

}

void ugrid_writer::nc_chk_ret(int status)
{
    if (status != NC_NOERR)
    {
        SPDLOG_ERROR("NC error status = {}", status);
        CHM_THROW_EXCEPTION(chm_error, nc_strerror(status));
    }
}

void ugrid_writer::close_ugrid()
{
    SPDLOG_DEBUG("Closing ugrid file");
    if (_ugrid_fid != -1)
    {
        nc_chk_ret(nc_close(_ugrid_fid));
    }
    _ugrid_fid = -1;

    // start writing at the start of new file
    _time_index = 0;

    // clean up for a new write
    _fname = "";
    _ugrid_id_var = {};

}
void ugrid_writer::write_ugrid(const std::vector<std::string>& output_variables)
{
    auto variables = output_variables.size() == 0 ? _mesh->face(0)->variables() : output_variables;
    if (_ugrid_fid==-1)
    {
        SPDLOG_DEBUG(_fname);
        // if we are resuming from checkpoint, don't mangle out existing ugrid!
        if (_global->from_checkpoint())
        {
            open_ugrid(variables);
        }
        else
        {
            init_ugrid(variables);
        }
    }

    // use C api as boost doesn't have info
    MPI_Comm comm = _comm_world;
    MPI_Info info_used;
    MPI_Comm_get_info(comm, &info_used);

    std::vector<size_t> all_face_offsets;
    boost::mpi::all_gather(_comm_world, _mesh->size_local_faces(), all_face_offsets);
    size_t offset_face = std::accumulate(all_face_offsets.begin(), all_face_offsets.begin() + static_cast<size_t>(_comm_world.rank()), 0);

    double time = _global->posix_time_double()  / 60 ; //s to m
    nc_chk_ret(nc_put_var1_double(_ugrid_fid, _ugrid_id_var["time"], &_time_index, &time));

    for (auto& var : variables)
    {
        std::vector<double> v(_mesh->size_local_faces());
        for (size_t i = 0; i < _mesh->size_local_faces(); i++)
        {
            double value = (*_mesh->face(i))[var];
            if (value == -9999.) value = nan("");
            v.at(i) = value;
        }

        size_t start[2] = {_time_index, offset_face};
        size_t count[2] = {1, _mesh->size_local_faces()};

        nc_chk_ret(nc_put_vara_double(_ugrid_fid, _ugrid_id_var[var], start,count,v.data()));

    }

    ++_time_index;

    MPI_Info_free(&info_used);
}
void ugrid_writer::open_ugrid(const std::vector<std::string>& output_variables)
{
    if (_ugrid_fid != -1)
    {
        CHM_THROW_EXCEPTION(model_init_error, "Netcdf ugrid file is already open");
    }

    if (!boost::filesystem::exists(_fname))
    {
        CHM_THROW_EXCEPTION(model_init_error, "Netcdf ugrid file is not found");
    }

    MPI_Comm comm = _comm_world;
    MPI_Info info_used;
    MPI_Comm_get_info(comm, &info_used);

    SPDLOG_DEBUG("Opening existing ugrid for writting");
    nc_chk_ret(nc_open_par(_fname.c_str(), NC_WRITE, _comm_world, info_used, &_ugrid_fid));

    for (const auto& vara:output_variables)
    {
        int status =
            nc_inq_varid(_ugrid_fid, vara.c_str(), &_ugrid_id_var[std::string(vara)]);
        if (status != NC_NOERR)
            CHM_THROW_EXCEPTION(model_init_error, "Netcdf ugrid file does not have variable to write: " + vara);
    }

    int status = nc_inq_varid(_ugrid_fid, "time", &_ugrid_id_var["time"]);
    if (status != NC_NOERR)
        CHM_THROW_EXCEPTION(model_init_error, "Netcdf ugrid file does not have variable to write: time");

    for (const auto& p : _ugrid_id_var)
    {
        nc_chk_ret(nc_var_par_access(_ugrid_fid, p.second, NC_COLLECTIVE));
    }

    int unlimdimidp;

    nc_chk_ret(nc_inq_unlimdim(_ugrid_fid, &unlimdimidp)); // get the time /dimension/. it's the only unlimited
    nc_chk_ret(nc_inq_dimlen(_ugrid_fid, unlimdimidp, &_time_index));

    SPDLOG_DEBUG("Existing ugrid output has {} timesteps already", _time_index);

    MPI_Info_free(&info_used);
}
void ugrid_writer::init_ugrid(const std::vector<std::string>& output_variables)
{
    timer c;
    // use C api as boost doesn't have info
    MPI_Comm comm = _comm_world;
    MPI_Info info_used;
    MPI_Comm_get_info(comm, &info_used);

    int status = nc_create_par(_fname.c_str(), NC_NETCDF4 | NC_CLOBBER, comm, info_used, &_ugrid_fid);
    if (status != NC_NOERR)
    {
        CHM_THROW_EXCEPTION(file_write_error, "Failed to create ugrid output file");
    }

    // We also only have per-rank vertex IDs (they aren't global)
    // thus need to compute a perrank offset to write into the global ugrid datastruct.
    // Rank 0 => [0, _num_local_vertex_on_rank0)
    // Rank 1 => [_num_local_vertex_on_rank0, _num_local_vertex_on_rank1)
    // etc

    SPDLOG_DEBUG("Starting ugrid def write");
    c.tic();

    std::vector<size_t> all_offsets;
    boost::mpi::all_gather(_comm_world, _mesh->size_local_vertex(), all_offsets);

    // Our rank needs the sum of all the vertexes from the preceeding ranks as its start offset
    size_t offset = std::accumulate(all_offsets.begin(), all_offsets.begin() + static_cast<size_t>(_comm_world.rank()), 0);

    std::vector<size_t> all_face_offsets;
    boost::mpi::all_gather(_comm_world, _mesh->size_local_faces(), all_face_offsets);
    size_t offset_face = std::accumulate(all_face_offsets.begin(), all_face_offsets.begin() + static_cast<size_t>(_comm_world.rank()), 0);

    // When the mesh is read we only know the number of vertexes on each rank, not the global number
    // This gathers the number of vertexes per rank.
    size_t num_global_vertex;
    boost::mpi::all_reduce(_comm_world, _mesh->size_local_vertex(), num_global_vertex, std::plus<size_t>());

    // Time-dependent face writes must align chunks with rank ownership to stay contiguous.
    size_t max_faces_per_rank = 0;
    boost::mpi::all_reduce(_comm_world, _mesh->size_local_faces(), max_faces_per_rank, boost::mpi::maximum<size_t>());
    if (_mesh->size_global_faces() > 0 && max_faces_per_rank > _mesh->size_global_faces())
        max_faces_per_rank = _mesh->size_global_faces();
    if (max_faces_per_rank == 0)
        max_faces_per_rank = 1;

    // Paraview's vtk-based ugrid reader segfaults when the mesh indexing is uint64. We probably don't have
    if (_mesh->size_global_faces() > UINT_MAX)
    {
        CHM_THROW_EXCEPTION(mesh_error, "Due to a limitation for paraview, the ugrid field can't have more than UINT_MAX nodes");
    }

    // base dimension structure
    int dim_Mesh2_node, dim_Mesh2_face, dim_two, dim_three;
    nc_chk_ret(nc_def_dim(_ugrid_fid, "nMesh2_node", num_global_vertex, &dim_Mesh2_node));
    nc_chk_ret(nc_def_dim(_ugrid_fid, "nMesh2_face", _mesh->size_global_faces(), &dim_Mesh2_face));
    nc_chk_ret(nc_def_dim(_ugrid_fid, "Two", 2, &dim_two));
    nc_chk_ret(nc_def_dim(_ugrid_fid, "Three", 3, &dim_three));

    // time Dimension
    int time_dimid, time_varid;

    // nc_def_dim(_ugrid_fid, "time", _global->n_timesteps(), &time_dimid);
    nc_chk_ret(nc_def_dim(_ugrid_fid, "time", NC_UNLIMITED, &time_dimid));
    nc_chk_ret(nc_def_var(_ugrid_fid, "time", NC_DOUBLE, 1, &time_dimid, &time_varid));
    size_t time_chunk[1] = {1};
    nc_chk_ret(nc_def_var_chunking(_ugrid_fid, time_varid, NC_CHUNKED, time_chunk));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, time_varid, "standard_name", strlen("time"), "time"));
    nc_chk_ret( nc_put_att_text(_ugrid_fid, time_varid, "long_name", strlen("Time"), "Time"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, time_varid, "units", strlen("minutes since 1970-01-01 00:00:00"), "minutes since 1970-01-01 00:00:00"));

    int var_Mesh2, var_Mesh2_face_nodes, var_Mesh2_node_x, var_Mesh2_node_y, var_Mesh2_node_z, var_Mesh2_node_z_PV;
    int dims_face_nodes[2] = {dim_Mesh2_face, dim_three};
    int dims_node[1] = {dim_Mesh2_node};

    // Mesh2 variable
    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2", NC_INT, 0, NULL, &var_Mesh2));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2, "cf_role", strlen("mesh_topology"), "mesh_topology"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2, "long_name", strlen("Topology data of 2D unstructured mesh"), "Topology data of 2D unstructured mesh"));

    int topo_dim = 2;
    nc_chk_ret(nc_put_att_int(_ugrid_fid, var_Mesh2, "topology_dimension", NC_INT, 1, &topo_dim));

    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2, "node_coordinates", strlen("Mesh2_node_x Mesh2_node_y"), "Mesh2_node_x Mesh2_node_y"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2, "face_node_connectivity", strlen("Mesh2_face_nodes"), "Mesh2_face_nodes"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2, "face_dimension", strlen("nMesh2_face"), "nMesh2_face"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2, "face_coordinates", strlen("Mesh2_face_x Mesh2_face_y"), "Mesh2_face_x Mesh2_face_y"));

    // Mesh2_face_nodes node connectivity that makes up the faces
    // should be NC_UINT64 but this crashes paraview
    // https://gitlab.kitware.com/paraview/paraview/-/issues/23019
    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_face_nodes", NC_UINT, 2, dims_face_nodes, &var_Mesh2_face_nodes));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_nodes, "cf_role", strlen("face_node_connectivity"), "face_node_connectivity"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_nodes, "long_name",
        strlen("Maps every triangular face to its three corner nodes."), "Maps every triangular face to its three corner nodes."));
    // UGRID recommends start_index attribute; indices are 0-based here
    {
        int start_index = 0;
        nc_chk_ret(nc_put_att_int(_ugrid_fid, var_Mesh2_face_nodes, "start_index", NC_INT, 1, &start_index));
    }

    // Mesh2_node_x
    double nan_value = NAN; // IEEE NaN
    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_node_x", NC_DOUBLE, 1, dims_node, &var_Mesh2_node_x));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_x, "standard_name", strlen("longitude"), "longitude"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_x, "long_name", strlen("Longitude of 2D mesh nodes."), "Longitude of 2D mesh nodes."));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_x, "units", strlen("degrees_east"), "degrees_east"));
    nc_chk_ret(nc_put_att_double(_ugrid_fid, var_Mesh2_node_x, "_FillValue", NC_DOUBLE, 1, &nan_value));

    // Mesh2_node_y
    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_node_y", NC_DOUBLE, 1, dims_node, &var_Mesh2_node_y));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_y, "standard_name", strlen("latitude"), "latitude"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_y, "long_name", strlen("Latitude of 2D mesh nodes."), "Latitude of 2D mesh nodes."));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_y, "units", strlen("degrees_north"), "degrees_north"));
    nc_chk_ret(nc_put_att_double(_ugrid_fid, var_Mesh2_node_y, "_FillValue", NC_DOUBLE, 1, &nan_value));

    // Mesh2_node_z elevation
    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_node_z", NC_DOUBLE, 1, dims_node, &var_Mesh2_node_z));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z, "standard_name", strlen("altitude"), "altitude"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z, "long_name", strlen("Z coordinate of 2D mesh nodes."), "Z coordinate of 2D mesh nodes."));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z, "units", strlen("m"), "m"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z, "mesh", strlen("Mesh2"), "Mesh2"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z, "location", strlen("node"), "node"));
    nc_chk_ret(nc_put_att_double(_ugrid_fid, var_Mesh2_node_z, "_FillValue", NC_DOUBLE, 1, &nan_value));

    // Scaled Mesh2_node_z elevation
    // Paraview struggles to plot the z coord when the x and y are in geographic, so this scales down the z
    // so it renders correctly.
    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_node_z_paraview", NC_DOUBLE, 1, dims_node, &var_Mesh2_node_z_PV));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z_PV, "standard_name", strlen("scaled_altitude"), "scaled_altitude"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z_PV, "long_name",
        strlen("Scaled Z coordinate of 2D mesh nodes."), "Scaled Z coordinate of 2D mesh nodes."));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z_PV, "units", strlen("m"), "m"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z_PV, "mesh", strlen("Mesh2"), "Mesh2"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_node_z_PV, "location", strlen("node"), "node"));
    nc_chk_ret(nc_put_att_double(_ugrid_fid, var_Mesh2_node_z_PV, "_FillValue", NC_DOUBLE, 1, &nan_value));


    int var_global_id, var_local_id, var_Mesh2_face_x, var_Mesh2_face_y, var_Mesh2_face_z;

    nc_chk_ret(nc_def_var(_ugrid_fid, "global_id", NC_UINT64, 1, &dim_Mesh2_face, &var_global_id));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_global_id, "mesh", strlen("Mesh2"), "Mesh2"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_global_id, "location", strlen("face"), "face"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_global_id, "coordinates", strlen("Mesh2_face_x Mesh2_face_y"), "Mesh2_face_x Mesh2_face_y"));

    nc_chk_ret(nc_def_var(_ugrid_fid, "local_id", NC_UINT64, 1, &dim_Mesh2_face, &var_local_id));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_local_id, "mesh", strlen("Mesh2"), "Mesh2"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_local_id, "location", strlen("face"), "face"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_local_id, "coordinates", strlen("Mesh2_face_x Mesh2_face_y"), "Mesh2_face_x Mesh2_face_y"));

    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_face_x", NC_DOUBLE, 1, &dim_Mesh2_face, &var_Mesh2_face_x));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_x, "standard_name", strlen("longitude"), "longitude"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_x, "long_name",
        strlen("Characteristics longitude of 2D mesh triangle (e.g. circumcenter coordinate)."),
        "Characteristics longitude of 2D mesh triangle (e.g. circumcenter coordinate)."));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_x, "units", strlen("degrees_east"), "degrees_east"));

    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_face_y", NC_DOUBLE, 1, &dim_Mesh2_face, &var_Mesh2_face_y));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_y, "standard_name", strlen("latitude"), "latitude"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_y, "long_name",
        strlen("Characteristics latitude of 2D mesh triangle (e.g. circumcenter coordinate)."),
        "Characteristics latitude of 2D mesh triangle (e.g. circumcenter coordinate)."));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_y, "units", strlen("degrees_north"), "degrees_north"));

    nc_chk_ret(nc_def_var(_ugrid_fid, "Mesh2_face_z", NC_DOUBLE, 1, &dim_Mesh2_face, &var_Mesh2_face_z));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_z, "standard_name", strlen("altitude"), "altitude"));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_z, "long_name",
        strlen("Characteristics latitude of 2D mesh triangle (e.g. circumcenter coordinate)."),
        "Characteristics latitude of 2D mesh triangle (e.g. circumcenter coordinate)."));
    nc_chk_ret(nc_put_att_text(_ugrid_fid, var_Mesh2_face_z, "units", strlen("m"), "m"));


    // bitgroom sigfigs
    int NSD = 4;

    int dims[2] = {time_dimid, dim_Mesh2_face};
    for (auto& var : output_variables)
    {
        nc_chk_ret(nc_def_var(_ugrid_fid, var.c_str(), NC_DOUBLE, 2, dims, &_ugrid_id_var[var]));
        size_t face_chunks[2] = {1, max_faces_per_rank};
        nc_chk_ret(nc_def_var_chunking(_ugrid_fid, _ugrid_id_var[var], NC_CHUNKED, face_chunks));
        nc_chk_ret(nc_put_att_text(_ugrid_fid, _ugrid_id_var[var], "mesh", strlen("Mesh2"), "Mesh2"));
        nc_chk_ret(nc_put_att_text(_ugrid_fid, _ugrid_id_var[var], "location", strlen("face"), "face"));
        nc_chk_ret(nc_put_att_double(_ugrid_fid, _ugrid_id_var[var], "_FillValue", NC_DOUBLE, 1, &nan_value));

        nc_chk_ret(nc_var_par_access(_ugrid_fid, _ugrid_id_var[var], NC_COLLECTIVE));

        if (compress) nc_chk_ret(nc_def_var_deflate(_ugrid_fid, _ugrid_id_var[var], 1, 1, 5));
        if (bitgroom) nc_chk_ret(nc_def_var_quantize(_ugrid_fid, _ugrid_id_var[var], NC_QUANTIZE_BITGROOM, NSD));
    }

    _ugrid_id_var["time"] = time_varid;
    nc_chk_ret(nc_var_par_access(_ugrid_fid, _ugrid_id_var["time"], NC_COLLECTIVE));
    if (compress) nc_chk_ret(nc_def_var_deflate(_ugrid_fid, _ugrid_id_var["time"], 1, 1, 5));


    std::map<std::string, int> param_id;
    if(_write_parameters)
    {
        auto define_face_variable = [&](const std::string& v, const std::string& prefix="") {
            return [&]() {
                nc_chk_ret(nc_def_var(_ugrid_fid, (prefix+v).c_str(), NC_DOUBLE, 1, &dim_Mesh2_face, &param_id[v]));
                size_t param_chunks[1] = {max_faces_per_rank};
                nc_chk_ret(nc_def_var_chunking(_ugrid_fid, param_id[v], NC_CHUNKED, param_chunks));
                nc_chk_ret(nc_put_att_text(_ugrid_fid, param_id[v], "mesh", strlen("Mesh2"), "Mesh2"));
                nc_chk_ret(nc_put_att_text(_ugrid_fid, param_id[v], "location", strlen("face"), "face"));
                nc_chk_ret(nc_put_att_double(_ugrid_fid, param_id[v], "_FillValue", NC_DOUBLE, 1, &nan_value));

                nc_chk_ret(nc_var_par_access(_ugrid_fid, param_id[v], NC_COLLECTIVE));

                if (compress) nc_chk_ret(nc_def_var_deflate(_ugrid_fid, param_id[v], 1, 1, 5));
                if (bitgroom) nc_chk_ret(nc_def_var_quantize(_ugrid_fid, param_id[v], NC_QUANTIZE_BITGROOM, NSD));
            };
        };


        auto params = _mesh->face(0)->parameters();
        for (auto &v: params)
        {
            define_face_variable(v, "param_")();
        }

        define_face_variable("Elevation")();
        define_face_variable("Slope")();
        define_face_variable("Aspect")();
        define_face_variable("Area")();


        nc_chk_ret(nc_def_var(_ugrid_fid, "owner", NC_INT, 1, &dim_Mesh2_face, &param_id["owner"]));
        size_t owner_chunks[1] = {max_faces_per_rank};
        nc_chk_ret(nc_def_var_chunking(_ugrid_fid, param_id["owner"], NC_CHUNKED, owner_chunks));
        nc_chk_ret(nc_put_att_text(_ugrid_fid, param_id["owner"], "mesh", strlen("Mesh2"), "Mesh2"));
        nc_chk_ret(nc_put_att_text(_ugrid_fid, param_id["owner"], "location", strlen("face"), "face"));
        nc_chk_ret(nc_var_par_access(_ugrid_fid, param_id["owner"], NC_COLLECTIVE));
        if (compress) nc_chk_ret(nc_def_var_deflate(_ugrid_fid, param_id["owner"], 1, 1, 5));

    }


    nc_chk_ret(nc_enddef(_ugrid_fid)); // End define mode

    auto t = c.toc<ms>();
    SPDLOG_DEBUG("Finished ugrid def section -- {} ms", t);


    // Fill node_x, node_y, node_z, face_nodes, face_neighbors, face_gid, static_param
    std::map<size_t, size_t> global_to_local_vertex_id;
    std::vector<size_t> global_vertex_id;

    OGRSpatialReference outsrs;
    outsrs.SetWellKnownGeogCS("WGS84");
    outsrs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER); //enforce ingoing as x y

    OGRSpatialReference insrs;
    insrs.importFromProj4(_mesh->proj4().c_str());
    insrs.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER); //enforce outgoing as x y

    OGRCoordinateTransformation* coordTrans = OGRCreateCoordinateTransformation(&insrs, &outsrs);

    c.tic();

    // write the vertexes
    {
        std::vector<double> v_x(_mesh->size_local_vertex());
        std::vector<double> v_y(_mesh->size_local_vertex());
        std::vector<double> v_z(_mesh->size_local_vertex());
        std::vector<double> v_z_scaled(_mesh->size_local_vertex());

        for (size_t i = 0; i < _mesh->size_local_vertex(); i++)
        {
            auto vit = _mesh->vertex(i);
            v_x.at(i) = vit->point().x();
            v_y.at(i) = vit->point().y();
            v_z.at(i) = vit->point().z();

            // scale Z for paraview 3D view
            v_z_scaled.at(i) = vit->point().z() / 100000. ;
        }

        if (!coordTrans->Transform(_mesh->size_local_vertex(), v_x.data(), v_y.data()))
        {
            CHM_THROW_EXCEPTION(forcing_error, "Failed to reproject coordinates");
        }

        size_t start_v[1] = {offset};
        size_t count_v[1] = {_mesh->size_local_vertex()};

        nc_chk_ret(nc_put_vara_double(_ugrid_fid, var_Mesh2_node_x, start_v, count_v, v_x.data()));
        nc_chk_ret(nc_put_vara_double(_ugrid_fid, var_Mesh2_node_y, start_v, count_v, v_y.data()));
        nc_chk_ret(nc_put_vara_double(_ugrid_fid, var_Mesh2_node_z, start_v, count_v, v_z.data()));

        nc_chk_ret(nc_put_vara_double(_ugrid_fid, var_Mesh2_node_z_PV, start_v, count_v, v_z_scaled.data()));
    }
    t = c.toc<ms>();
    SPDLOG_DEBUG("Finished ugrid vertex -- {} ms", t);


    c.tic();
    // Write the connectivity matrix that defines what vertexes each face is comprised of
    {
        boost::multi_array<unsigned int,2> connectivity(boost::extents[_mesh->size_local_faces()][3]);
        for (size_t i = 0; i < _mesh->size_local_faces(); i++)
        {
            auto fit = _mesh->face(i);
            connectivity[i][0] = static_cast<unsigned int>(fit->vertex(0)->get_id() + offset);
            connectivity[i][1] = static_cast<unsigned int>(fit->vertex(1)->get_id() + offset);
            connectivity[i][2] = static_cast<unsigned int>(fit->vertex(2)->get_id() + offset);
        }

        size_t start[2] = {offset_face, 0};   // Start at face_idx, first node
        size_t count[2] = {_mesh->size_local_faces(), 3};          // One face, all three nodes
        nc_chk_ret(nc_put_vara_uint(_ugrid_fid, var_Mesh2_face_nodes, start, count, connectivity.data()));
    }
    t = c.toc<ms>();
    SPDLOG_DEBUG("Finished ugrid connectivity -- {} ms", t);

    c.tic();
    // Write the x,y,z for the face centers
    {
        std::vector<double> f_x(_mesh->size_local_faces());
        std::vector<double> f_y(_mesh->size_local_faces());
        std::vector<double> f_z(_mesh->size_local_faces());

        for (size_t i = 0; i < _mesh->size_local_faces(); i++)
        {
            auto fit = _mesh->face(i);

           f_x.at(i) = fit->center().x();
           f_y.at(i) = fit->center().y();
           f_z.at(i) = fit->center().z();
        }

        if (!coordTrans->Transform(_mesh->size_local_faces(), f_x.data(), f_y.data()))
        {
            CHM_THROW_EXCEPTION(forcing_error, "Failed to reproject coordinates");
        }

        size_t start[1] = {offset_face};
        size_t count[1] = {_mesh->size_local_faces()};

        nc_chk_ret(nc_put_vara_double(_ugrid_fid, var_Mesh2_face_x, start, count, f_x.data()));
        nc_chk_ret(nc_put_vara_double(_ugrid_fid, var_Mesh2_face_y, start, count, f_y.data()));
        nc_chk_ret(nc_put_vara_double(_ugrid_fid, var_Mesh2_face_z, start, count, f_z.data()));
    }
    t = c.toc<ms>();
    SPDLOG_DEBUG("Finished ugrid xyz centres -- {} ms", t);

    c.tic();
    {
        size_t start[1] = {offset_face};
        size_t count[1] = {_mesh->size_local_faces()};

        // update to uing / size_t nc_put_var1_ulonglong
        nc_put_vara_int(_ugrid_fid, var_global_id, start, count, _mesh->get_global_IDs().data());
    }
    t = c.toc<ms>();
    SPDLOG_DEBUG("Finished ugrid global_id -- {} ms", t);

    c.tic();
    if(_write_parameters)
    {
        // This strategy takes a bit more CPU because of param* n_face iterations, but cuts down on having to store
        // entire duplicates of the mesh
        for (auto &v: _mesh->face(0)->parameters())
        {
            std::vector<double> param(_mesh->size_local_faces());

            for (size_t i = 0; i < _mesh->size_local_faces(); i++)
            {
                double p = _mesh->face(i)->parameter(v);
                if( p == -9999.) p = nan("");
                param.at(i) = p;
            }

            size_t start[1] = {offset_face};
            size_t count[1] = {_mesh->size_local_faces()};
            nc_chk_ret(nc_put_vara_double(_ugrid_fid, param_id[v], start, count, param.data()));
        }

        // Slope
        {
            std::vector<double> param(_mesh->size_local_faces());

            for (size_t i = 0; i < _mesh->size_local_faces(); i++)
            {
                double p = _mesh->face(i)->slope();
                if( p == -9999.) p = nan("");
                param.at(i) = p;
            }

            size_t start[1] = {offset_face};
            size_t count[1] = {_mesh->size_local_faces()};
            nc_chk_ret(nc_put_vara_double(_ugrid_fid, param_id["Slope"], start, count, param.data()));
        }

        // Aspect
        {
            std::vector<double> param(_mesh->size_local_faces());

            for (size_t i = 0; i < _mesh->size_local_faces(); i++)
            {
                double p = _mesh->face(i)->aspect();
                if( p == -9999.) p = nan("");
                param.at(i) = p;
            }

            size_t start[1] = {offset_face};
            size_t count[1] = {_mesh->size_local_faces()};
            nc_chk_ret(nc_put_vara_double(_ugrid_fid, param_id["Aspect"], start, count, param.data()));
        }

        // Area
        {
            std::vector<double> param(_mesh->size_local_faces());

            for (size_t i = 0; i < _mesh->size_local_faces(); i++)
            {
                double p = _mesh->face(i)->get_area();
                if( p == -9999.) p = nan("");
                param.at(i) = p;
            }

            size_t start[1] = {offset_face};
            size_t count[1] = {_mesh->size_local_faces()};
            nc_chk_ret(nc_put_vara_double(_ugrid_fid, param_id["Area"], start, count, param.data()));
        }


        // Elevation
        {
            std::vector<double> param(_mesh->size_local_faces());

            for (size_t i = 0; i < _mesh->size_local_faces(); i++)
            {
                double p = _mesh->face(i)->get_z();
                if( p == -9999.) p = nan("");
                param.at(i) = p;
            }

            size_t start[1] = {offset_face};
            size_t count[1] = {_mesh->size_local_faces()};
            nc_chk_ret(nc_put_vara_double(_ugrid_fid, param_id["Elevation"], start, count, param.data()));
        }

        // Slope
        {
            std::vector<int> param(_mesh->size_local_faces(), _comm_world.rank());

            size_t start[1] = {offset_face};
            size_t count[1] = {_mesh->size_local_faces()};
            nc_chk_ret(nc_put_vara_int(_ugrid_fid, param_id["owner"], start, count, param.data()));
        }
    }
    t = c.toc<ms>();
    SPDLOG_DEBUG("Finished ugrid params -- {} ms", t);

    // clean up the handle to info
    MPI_Info_free(&info_used);

    OGRCoordinateTransformation::DestroyCT(coordTrans);

}
