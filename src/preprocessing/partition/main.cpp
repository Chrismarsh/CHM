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

#include "H5Cpp.h"
#include "core.hpp"
#include "logger.hpp"
#include "sort_perm.hpp"
#include "triangulation.hpp"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>

#include <algorithm>
#include <numeric>
#include <string>

namespace pt = boost::property_tree;
namespace po = boost::program_options;

using namespace H5;


#define USE_MPI 1

typedef struct _MeshParameters
{
    std::vector<std::string> names;
} MeshParameters;

// this is declared in the main triangulation.cpp
herr_t group_info(hid_t loc_id, const char* name, const H5L_info_t* linfo, void* opdata);

class preprocessingTriangulation : public triangulation
{
  public:
    preprocessingTriangulation()
    {
        _is_standalone = false;
        _is_partition = true;
        _write_ghost_neighbors_to_vtu = true;
        _write_parameters_to_vtu = false;
    }

    void write_vtu(std::string file_name, std::vector<std::string> output_variables ={} )
    {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();

        triangles->Allocate(_local_faces.size());

        vtkSmartPointer<vtkStringArray> proj4 = vtkSmartPointer<vtkStringArray>::New();
        proj4->SetNumberOfComponents(1);
        proj4->SetName("proj4");
        proj4->InsertNextValue(_srs_wkt);

        double scale = is_geographic() == true ? 100000. : 1.;

        std::map<int, int> global_to_local_vertex_id;
        std::vector<int> global_vertex_id;

        // npoints holds the total number of points
        int npoints=0;
        for (size_t i = 0; i < _local_faces.size(); i++)
        {
            mesh_elem fit = _local_faces.at(i);

            vtkSmartPointer<vtkTriangle> tri =
                vtkSmartPointer<vtkTriangle>::New();

            // loop over vertices of a face
            for (int j=0;j<3;++j)
            {
                auto vit = fit->vertex(j);
                int global_id = vit->get_id();
                // If point hasn't been seen yet, account for it
                if ( global_to_local_vertex_id.find(global_id) == global_to_local_vertex_id.end() ) {
                    global_to_local_vertex_id[global_id] = npoints;
                    npoints++;
                    points->InsertNextPoint(vit->point().x()*scale, vit->point().y()*scale, vit->point().z());
                    global_vertex_id.push_back(global_id);
                }
                tri->GetPointIds()->SetId(j, global_to_local_vertex_id[global_id]);
            }

            triangles->InsertNextCell(tri);
        }
        _vtk_unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        _vtk_unstructuredGrid->SetPoints(points);
        _vtk_unstructuredGrid->SetCells(VTK_TRIANGLE, triangles);
        _vtk_unstructuredGrid->GetFieldData()->AddArray(proj4);

        auto variables = output_variables.size() == 0 ? this->face(0)->variables() : output_variables;
        for(auto& v: variables)
        {
            data[v] = vtkSmartPointer<vtkFloatArray>::New();
            data[v]->SetName(v.c_str());
        }

        data["Elevation"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Elevation"]->SetName("Elevation");

        // Global vertex ids -> only need to be set here, get written in the writer
        vertex_data["global_id"] = vtkSmartPointer<vtkFloatArray>::New();
        vertex_data["global_id"]->SetName("global_id");
        for(int i=0;i<npoints;++i){
            vertex_data["global_id"]->InsertTuple1(i,global_vertex_id[i]);
        }

        for (size_t i = 0; i < _local_faces.size(); i++)
        {
            mesh_elem fit = _local_faces.at(i);

            for (auto &v: variables)
            {
                double d = (*fit)[v];
                if(d == -9999.)
                {
                    d = nan("");
                }

                data[v]->InsertTuple1(i,d);
            }

            data["Elevation"]->InsertTuple1(i,fit->get_z());

        }
        for(auto& m : data)
        {
            _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);
        }

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(file_name.c_str());

        writer->SetInputData(_vtk_unstructuredGrid);

        writer->Write();
    }

    // sets the partition owner on every face
    // we don't do a valid rank check here as every triangle will need it set so ghosts are correctly ID'd
    void set_face_partition_owner(int mpi_rank)
    { // Set up so that all processors know how 'big' all other processors are
        if (_comm_world.size() == 1)
        {
            _num_faces_in_partition.resize(1);
            _num_faces_in_partition.at(0) = _num_global_faces;
        }
        else
        {
            _num_faces_in_partition.resize(_comm_world.size());
            for (size_t i = 0; i < _comm_world.size(); ++i)
            {
                _num_faces_in_partition.at(i) = _local_sizes.at(i);
            }
        }
        // each processor only knows its own start and end indices
        size_t face_start_idx = 0;
        size_t face_end_idx = _num_faces_in_partition.at(0) - 1;
        for (int i = 1; i <= _comm_world.rank(); ++i)
        {
            face_start_idx += _num_faces_in_partition.at(i - 1);
            face_end_idx += _num_faces_in_partition.at(i);
        }

#pragma omp parallel for
        for (size_t local_ind = 0; local_ind < _num_faces_in_partition.at(_comm_world.rank()); ++local_ind)
        {
            size_t offset_idx = face_start_idx + local_ind;
            _faces.at(_global_IDs.at(offset_idx))->owner = mpi_rank;
        }

    }

    /**
     * This partitions the domain for the given rank. Differs from `partition_mesh` in triangulation by:
     * - don't set the face->global_id as that is already set by set_face_parition_owner
     * - don't set the global_IDs vector as that is read from the input h5 files
     * @param mpi_rank
     */
    void partition(int mpi_rank)
    { // Set up so that all processors know how 'big' all other processors are
        if (_comm_world.size() == 1)
        {
            _num_faces_in_partition.resize(1);
            _num_faces_in_partition.at(0) = _num_global_faces;
        }
        else
        {
            _num_faces_in_partition.resize(_comm_world.size());
            for (size_t i = 0; i < _comm_world.size(); ++i)
            {
                _num_faces_in_partition.at(i) = _local_sizes.at(i);
            }
        }
        // each processor only knows its own start and end indices
        size_t face_start_idx = 0;
        size_t face_end_idx = _num_faces_in_partition.at(0) - 1;
        for (int i = 1; i <= _comm_world.rank(); ++i)
        {
            face_start_idx += _num_faces_in_partition.at(i - 1);
            face_end_idx += _num_faces_in_partition.at(i);
        }
        // Store in the public members
        global_cell_start_idx = face_start_idx;
        global_cell_end_idx = face_end_idx;

        // Set size of vector containing locally owned faces
        _local_faces.resize(_num_faces_in_partition.at(_comm_world.rank()));

        // Loop can't be parallel due to modifying map
        for (size_t local_ind = 0; local_ind < _local_faces.size(); ++local_ind)
        {
            size_t offset_idx = face_start_idx + local_ind;
            _global_to_locally_owned_index_map[_global_IDs.at(offset_idx)] = local_ind;

            _faces.at(_global_IDs.at(offset_idx))->is_ghost = false;
            _faces.at(_global_IDs.at(offset_idx))->cell_local_id = local_ind;
            _local_faces.at(local_ind) = _faces.at(_global_IDs.at(offset_idx));
        }

        LOG_DEBUG << "MPI Process " << mpi_rank << ": start " << face_start_idx << ", end " << face_end_idx
                  << ", number " << _local_faces.size();
    }
    void read_h5(const std::string& mesh_filename)
    {
        std::vector<Point_2> center_points;

        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            Exception::dontPrint();

            // Open an existing file and dataset.
            H5File file(mesh_filename, H5F_ACC_RDONLY);

            // check the mesh version first
            {
                std::string v;
                try
                {
                    // Read the mesh version
                    H5::DataSpace dataspace(1, &meshver_dims);
                    H5::Attribute attribute = file.openAttribute("/mesh/version");
                    attribute.read(meshver_t, v);
                }
                catch (AttributeIException& e)
                {
                    v = ""; // will cause a default to version 1.0.0
                }
                _version.from_string(v);

                if (!_version.mesh_ver_meets_min_partition())
                    CHM_THROW_EXCEPTION(mesh_error, "h5 mesh doesn't meet criteria to partition -- version to too old "
                                                    "or was not permuted with -t metis?");
            }

            std::vector<std::array<double, 3>> vertex;
            std::vector<std::array<int, 3>> elem;
            std::vector<std::array<int, 3>> neigh;

            {
                // Read the proj4
                H5::DataSpace dataspace(1, &proj4_dims);
                H5::Attribute attribute = file.openAttribute("/mesh/proj4");
                attribute.read(proj4_t, _srs_wkt);
            }

            { // Write the is_geographic
                H5::DataSpace dataspace(1, &geographic_dims);
                H5::Attribute attribute = file.openAttribute("/mesh/is_geographic");
                attribute.read(PredType::NATIVE_HBOOL, &_is_geographic);
            }

            try
            {
                {
                    // Read the proj4
                    H5::DataSpace dataspace(1, &partition_type_dims);
                    H5::Attribute attribute = file.openAttribute("/mesh/partition_method");
                    attribute.read(partition_type_t, _partition_method);
                }

                {
                    H5::DataSet dataset = file.openDataSet("/mesh/local_sizes");
                    H5::DataSpace dataspace = dataset.getSpace();
                    hsize_t nelem;
                    int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);
                    _local_sizes.resize(nelem);
                    dataset.read(_local_sizes.data(), PredType::NATIVE_INT);
                }
            }
            catch (...)
            {
                CHM_THROW_EXCEPTION(mesh_error,
                                    "This h5 was produced by an older version of partition/meshpermutation.py and is "
                                    "lacking a key field. Please rerun these tools");
            }

            {
                DataSet dataset = file.openDataSet("/mesh/cell_global_id");
                DataSpace dataspace = dataset.getSpace();
                hsize_t nelem;
                int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);
                _global_IDs.resize(nelem);
                dataset.read(_global_IDs.data(), PredType::NATIVE_INT);
            }

            {
                // Open the vertices dataset
                DataSet dataset = file.openDataSet("/mesh/vertex");
                DataSpace dataspace = dataset.getSpace();

                hsize_t nvert;
                int ndims = dataspace.getSimpleExtentDims(&nvert, NULL);

                // Ensure enough space in the vector
                vertex.resize(nvert);
                // Default args read all of the dataspace
                dataset.read(vertex.data(), vertex_t);

                for (size_t i = 0; i < nvert; i++)
                {

                    Point_3 pt(vertex[i][0], vertex[i][1], vertex[i][2]);
                    _max_z = std::max(_max_z, vertex[i][2]);
                    _min_z = std::min(_min_z, vertex[i][2]);

                    Vertex_handle Vh = create_vertex();
                    Vh->set_point(pt);
                    Vh->set_id(i);
                    _vertexes.push_back(Vh);
                }

                _num_vertex = _vertexes.size();
                LOG_DEBUG << "# nodes created = " << _num_vertex;
            }

            {
                DataSet dataset = file.openDataSet("/mesh/elem");
                DataSpace dataspace = dataset.getSpace();

                hsize_t nelem;
                int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);

                // Ensure enough space in the vector
                elem.resize(nelem);
                // Default args read all of the dataspace
                dataset.read(elem.data(), elem_t);

                for (size_t i = 0; i < nelem; i++)
                {

                    auto vert1 = _vertexes.at(elem[i][0]); // 0 indexing
                    auto vert2 = _vertexes.at(elem[i][1]);
                    auto vert3 = _vertexes.at(elem[i][2]);

                    auto face = create_face(vert1, vert2, vert3);
                    face->cell_global_id = i;
                    face->cell_local_id = i;

                    if (_is_geographic)
                    {
                        face->_is_geographic = true;
                    }

                    face->_debug_ID = -(i + 1); // all ids will be negative starting at -1. Named ids (for output) will
                                                // be positive starting at 0
                    face->_debug_name = std::to_string(i);

                    vert1->set_face(face);
                    vert2->set_face(face);
                    vert3->set_face(face);

                    _faces.push_back(face);
                }

                _num_faces = _faces.size();
                LOG_DEBUG << "Created a mesh with " << size_faces() << " triangles";
            }

            {

                DataSet dataset = file.openDataSet("/mesh/neighbor");
                DataSpace dataspace = dataset.getSpace();

                hsize_t nelem;
                int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);

                if (elem.size() != nelem)
                {
                    BOOST_THROW_EXCEPTION(config_error()
                                          << errstr_info("Expected: " + std::to_string(elem.size()) +
                                                         " neighborlists, got: " + std::to_string(nelem)));
                }

                // Ensure enough space in the vector
                neigh.resize(nelem);
                // Default args read all of the dataspace
                dataset.read(neigh.data(), neighbor_t);

                LOG_DEBUG << "Building face neighbors";
                for (size_t i = 0; i < nelem; i++)
                {
                    auto face = _faces.at(i);

                    if (neigh[i][0] > static_cast<int>(nelem) || neigh[i][1] > static_cast<int>(nelem) ||
                        neigh[i][2] > static_cast<int>(nelem))
                    {
                        BOOST_THROW_EXCEPTION(config_error() << errstr_info("Face " + std::to_string(i) +
                                                                            " has out of bound neighbors."));
                    }

                    //-1 is now the no neighbor value
                    Face_handle face0 = neigh[i][0] != -1 ? _faces.at(neigh[i][0]) : nullptr; // Face_handle()
                    Face_handle face1 = neigh[i][1] != -1 ? _faces.at(neigh[i][1]) : nullptr;
                    Face_handle face2 = neigh[i][2] != -1 ? _faces.at(neigh[i][2]) : nullptr;

                    face->set_neighbors(face0, face1, face2);
                }
            }

        } // end of try block

        // catch failure caused by the H5File operations
        catch (FileIException& error)
        {
            error.printErrorStack();
        }

        // catch failure caused by the DataSet operations
        catch (DataSetIException& error)
        {
            error.printErrorStack();
        }
    }

    /**
     * This differences from triangulation::determine_process_ghost_faces_nearest_neighbors by including a modification
     * to set the neigh ghost type
     */
    void determine_process_ghost_faces_nearest_neighbors()
    {
        // NOTE that this algorithm is not implemented for multithread
        // - multithread can be implemented similarly to determine_local_boundary_faces
        //   - not a priority since number of boundary faces should be small

        // Ensure that the local boundary faces have been determined, but ghost nearest neighbors have not been set
        assert(_boundary_faces.size() != 0);
        assert(_ghost_neighbors.size() == 0);

        LOG_DEBUG << "Determining ghost region info";

        // Vector for append speed
        std::vector<mesh_elem> ghosted_boundary_nearest_neighbors;

        for (size_t face_index = 0; face_index < _boundary_faces.size(); ++face_index)
        {
            // face_index is a local index... get the face handle
            auto face = _boundary_faces.at(face_index).first;
            // append the ghosted nearest neighbors
            for (int i = 0; i < 3; ++i)
            {
                auto neigh = face->neighbor(i);
                if (neigh != nullptr && neigh->is_ghost)
                {
                    ghosted_boundary_nearest_neighbors.push_back(neigh);
                    neigh->get_module_data<ghost_info>("partition_tool").ghost_type = ghost_info::GHOST_TYPE::NEIGH;
                }
            }
        }

        // Convert to a set to remove duplicates
        std::unordered_set<mesh_elem> tmp_set(std::begin(ghosted_boundary_nearest_neighbors),
                                              std::end(ghosted_boundary_nearest_neighbors));
        // Convert the set to a vector
        _ghost_neighbors.insert(std::end(_ghost_neighbors), std::begin(tmp_set), std::end(tmp_set));
        // Sort the ghost neighbors
        // NOTE:
        // - sorting this vector by cell_global_id effectively partitions the ghost neighbors to be contiguous in
        // communication partners
        // This is required to put them into chunked sections
        std::sort(_ghost_neighbors.begin(), _ghost_neighbors.end(),
                  [&](const auto& a, const auto& b) { return a->cell_global_id < b->cell_global_id; });

        LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _ghost_neighbors.size()
                  << " ghosted nearest neighbors.";
    }

    void from_file_and_partition(const std::string& mesh_filename, const std::vector<std::string>& param_filenames,
                                 size_t MPI_ranks, double max_ghost_distance, int standalone_rank,
                                 std::vector<int> ranks_to_keep, bool output_vtu)
    {

        _comm_world._size = MPI_ranks;

        {
            core c; // instantiate this long enough to use this one function
            _is_geographic = c.check_is_geographic(mesh_filename);
        }

        if (_is_geographic)
        {
            math::gis::point_from_bearing = &math::gis::point_from_bearing_latlong;
            math::gis::distance = &math::gis::distance_latlong;
        }
        else
        {
            math::gis::point_from_bearing = &math::gis::point_from_bearing_UTM;
            math::gis::distance = &math::gis::distance_UTM;
        }

        auto mesh_file_extension = boost::filesystem::path(mesh_filename).extension().string();

        if (mesh_file_extension == ".partition")
        {
            CHM_THROW_EXCEPTION(mesh_error, "Cannot run this tool with a .partition mesh input!");
        }

        if (mesh_file_extension == ".h5")
        {
            LOG_DEBUG << "Partitioning mesh " << mesh_filename << " with #ranks=" << MPI_ranks;
            read_h5(mesh_filename);
        }
        else
        {
            auto mesh = read_json(mesh_filename);

            // mesh will have been loaded by the geographic check so don't re load it here
            bool triarea_found = false;

            // Parameter files
            for (auto param_file : param_filenames)
            {
                pt::ptree param_json = read_json(param_file);

                for (auto& ktr : param_json)
                {
                    // use put to ensure there are no duplciate parameters...
                    std::string key = ktr.first.data();
                    mesh.put_child("parameters." + key, ktr.second);

                    if (key == "area")
                        triarea_found = true;

                    LOG_DEBUG << "Inserted parameter " << ktr.first.data() << " into the config tree.";
                }
            }

            if (_is_geographic && !triarea_found)
            {
                BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Geographic meshes require the triangle area be "
                                                                  "present in a .param file. Please include this."));
            }

            // Initial condition files
            //            for(auto ic_file : initial_condition_file_paths)
            //            {
            //                pt::ptree ic_json = read_json(ic_file);
            //
            //                for(auto& ktr : ic_json)
            //                {
            //                    //use put to ensure there are no duplciate parameters...
            //                    std::string key = ktr.first.data();
            //                    mesh.put_child( "initial_conditions." + key ,ktr.second);
            //                    LOG_DEBUG << "Inserted initial condition " << ktr.first.data() << " into the config
            //                    tree.";
            //                }
            //            }

            try
            {
                std::vector<size_t> permutation;
                mesh.get_child("mesh.cell_global_id");
            }
            catch (pt::ptree_bad_path& e)
            {

                LOG_WARNING
                    << "No face permutation was found in the mesh file. This will result in poor MPI performance"
                       " due to increased communications. \nPlease see the mesh permutation documentation for more "
                       "details\n"
                       "https://mesher-hydro.readthedocs.io/en/latest/tools.html#mesherpermuation-py";
            }
            from_json(mesh);

            LOG_DEBUG << "Converting mesh to hdf5...";
            auto base_name = boost::filesystem::path(mesh_filename).stem().string();
            triangulation::to_hdf5(base_name);
            LOG_DEBUG << "HDF5 written, terminating";
            return;
        }

        pt::ptree tree;
        tree.put("ranks", ranks_to_keep.size() == 0? MPI_ranks : ranks_to_keep.size());
        tree.put("max_ghost_distance", max_ghost_distance);

        pt::ptree meshes;
        pt::ptree params;

        std::vector<size_t> tmp_local_sizes;
        pt::ptree local_sizes;
        for (size_t i = 0; i < _local_sizes.size(); ++i)
        {
            bool valid_rank = (std::find(ranks_to_keep.begin(), ranks_to_keep.end(),i) != ranks_to_keep.end());
            if(!valid_rank)
                continue;

            pt::ptree s;
            tmp_local_sizes.push_back(_local_sizes.at(i));
            s.put("", std::to_string(_local_sizes.at(i)));
            local_sizes.push_back(std::make_pair("", s));
        }
        tree.add_child("local_sizes", local_sizes);

        // If we are not running a specific set of ranks, then this will be the same as _faces.size();
        _num_global_faces = std::accumulate(tmp_local_sizes.begin(), tmp_local_sizes.end(),
                                       decltype(tmp_local_sizes)::value_type(0));
        tree.put("num_global_faces", _num_global_faces);

        std::string filename_base = mesh_filename.substr(0, mesh_filename.length() - 3);

        auto partition_dir =
            boost::filesystem::path(filename_base + ".np" + std::to_string(_comm_world.size()) + ".partition.meshes");
        boost::filesystem::create_directory(partition_dir);

        int start_rank = 0;
        int end_rank = MPI_ranks;

        // if we output a standalone mesh only process that single mpirank
        if (standalone_rank > -1)
        {
            start_rank = standalone_rank;
            end_rank = standalone_rank + 1;
            _is_standalone = true;
        }

        if(output_vtu)
        {
            // init the datastructs to hold information for outputting to VTU
            std::set< std::string > vtu_outputs = { "owner", "is_ghost", "ghost_type", "global_id","local_id"};
#pragma omp parallel for
            for (size_t i = 0; i < _faces.size(); ++i)
            {
                auto f = _faces.at(i);
                f->init_time_series(vtu_outputs);
            }
        }


        // if we are only outputting a subset of the domain, we need to rewrite the owner to account for the number
        // of mpi ranks that will participate in the subset. E.g., if we have 500 ranks but only output 3:
        // -v 400 -v 100 -v 405
        // then 100 will be mapped to rank 0
        // 400 mapped to rank 1
        // 405 mapped to rank 2
        // for a total of 3 ranks for the 3 subsets
        size_t proxy_global_id = 0;
        std::map<size_t,size_t> _global_proxy_to_global_index_map;
        std::map<size_t,size_t> _global_id_to_proxy_index_map;

        bool is_subset = (MPI_ranks != ranks_to_keep.size());
        if(is_subset)
            LOG_WARNING << "Processing will produce only subset output because tool was invoked with -v flag(s)";

        // set all our owners first
        LOG_DEBUG << "Determining partition owners";

        std::map<size_t,size_t> _mpirank_proxy_to_mpirank_map;
        std::map<size_t,size_t> _mpirank_to_proxy_mpirank_map;
        int proxy_mpirank = start_rank;
        for (int mpirank = start_rank; mpirank < end_rank; mpirank++)
        {
            _comm_world._rank = mpirank;
            set_face_partition_owner(mpirank);

            bool valid_rank = (std::find(ranks_to_keep.begin(), ranks_to_keep.end(),mpirank) != ranks_to_keep.end());
            if(valid_rank)
            {
                _mpirank_proxy_to_mpirank_map[proxy_mpirank] = mpirank;
                _mpirank_to_proxy_mpirank_map[mpirank] = proxy_mpirank;

                LOG_DEBUG<< "MPI rank = " << mpirank << " will be remapped to  " << _mpirank_to_proxy_mpirank_map[mpirank];

                proxy_mpirank++;
            }
        }




        // iterate over all the ranks but since we might be only doing a subset, we need to maintain a proxy rank
        for (int mpirank = start_rank; mpirank < end_rank; mpirank++)
        {
            bool valid_rank = (std::find(ranks_to_keep.begin(), ranks_to_keep.end(),mpirank) != ranks_to_keep.end());
            if(!valid_rank)
                continue;

            _comm_world._rank = mpirank;

            // reset everything
            _global_to_locally_owned_index_map.clear();
            _num_faces_in_partition = {};
            _local_faces = {};
            _boundary_faces = {};
            _ghost_neighbors = {};
            _ghost_neighbor_owners = {};
            _ghost_faces = {};
            _comm_partner_ownership.clear();
            _global_index_to_local_ghost_map.clear();
            global_indices_to_send.clear();
            _param_data.clear();

            std::set<std::string> pt = {"partition_tool"};

#pragma omp parallel for
            for (size_t i = 0; i < _faces.size(); ++i)
            {
                auto f = _faces.at(i);
                f->is_ghost = true; // default state, switches to false later
                f->init_module_data(pt);

                auto& gi = f->make_module_data<ghost_info>("partition_tool");

                // doesn't match the is_ghost default state. Here we assume false, and then switch it to the correct
                // type when determined
                gi.ghost_type = ghost_info::GHOST_TYPE::NONE;

                if(output_vtu)
                {
                    (*f)["owner"] = -9999;
                    (*f)["is_ghost"] = -9999;
                    (*f)["ghost_type"]=-9999;
                    (*f)["global_id"]=-9999;
                    (*f)["local_id"]=-9999;
                }

            }

            partition(_comm_world._rank);

            _num_faces = _local_faces.size();

            if( _num_faces != _local_sizes.at(mpirank))
            {
                LOG_ERROR<< "Read in #faces=" << _num_faces <<" faces but the partition defined #faces=" << _local_sizes.at(mpirank);
                CHM_THROW_EXCEPTION(mesh_error, "Number of local faces does not match the expected number of faces defined in the partition");
            }

            // determine what faces in our partition are the boundary faces
            determine_local_boundary_faces();

            // figure out which neighbours of the  boundary faces are ghosts -- Type I ghosts
            determine_process_ghost_faces_nearest_neighbors();

            // TODO: Need to auto-determine how far to look based on module setups
            // when this is called, it will output
            // MPI Process 0 has XXX ghosted faces.
            // regardless of what MPI rank we are here as it is using the super's _commworld for the output /only/

            // Find type II ghosts
            determine_process_ghost_faces_by_distance(max_ghost_distance);

            // _ghost_faces as set by determine_process_ghost_faces_by_distance is missing the nearest neighbours
            // and so add them in here
            _ghost_faces.insert(std::end(_ghost_faces),
                                std::begin(_ghost_neighbors),std::end(_ghost_neighbors));

            // Convert to a set to remove duplicates
            std::unordered_set<mesh_elem> tmp_set(std::begin(_ghost_faces),
                                                  std::end(_ghost_faces));

            // rebuild after removing duplicates
            _ghost_faces.clear();
            _ghost_faces.insert(std::end(_ghost_faces),
                                std::begin(tmp_set),std::end(tmp_set));

            // sort as this is needed for the chunking of ghosts per mpi owner
            std::sort(_ghost_faces.begin(), _ghost_faces.end(),
                      [&](const auto& a, const auto& b) { return a->cell_global_id < b->cell_global_id; });

            // determine the MPI rank that owns our Type I ghosts
            determine_ghost_owners();


            // flag all the ghosts that aren't Type I (neigh) as Typde II distance ghosts
#pragma omp parallel for
            for (size_t i = 0; i < _ghost_faces.size(); i++)
            {
                auto face = _ghost_faces.at(i);
                auto& gi = face->get_module_data<ghost_info>("partition_tool");

                // if we aren't a neigh, we were added by determine_process_ghost_faces_by_distance
                if (gi.ghost_type != ghost_info::GHOST_TYPE::NEIGH)
                    gi.ghost_type = ghost_info::GHOST_TYPE::DIST;
            }


            // Determine if there are any ghost faces owned by an MPI rank that we are not outputting. Default case
            // is that there will be no invalid ghosts unless we use the -v option
            auto invalid_ghosts = std::remove_if(
                _ghost_faces.begin(), _ghost_faces.end(),
                [&](auto& itr)
                {
                    bool valid_rank =
                        (std::find(ranks_to_keep.begin(), ranks_to_keep.end(), itr->owner) != ranks_to_keep.end());
                    return !valid_rank;
                });
            _ghost_faces.erase(invalid_ghosts, _ghost_faces.end());

            // load params
            std::unique_ptr<MeshParameters> pars = read_h5_params(param_filenames);

            if (!_is_standalone)
            {
                // now, augment the local faces with our ghosts, so that the ghosts are
                // a) flagged and b) available in this subset
                // _ghost_faces includes /all/ the ghost faces: neighbour + distance
                _local_faces.insert(_local_faces.end(), _ghost_faces.begin(), _ghost_faces.end());
            }

            // if we are outputting a subset, shim in all the proxy ranks and global IDs here
            if(is_subset)
            {
                for(auto& f: _local_faces)
                {
                    // If we haevn't seen this ID, add it to our list and inc our global id counter
                    if(_global_id_to_proxy_index_map.find( f->cell_global_id) ==
                        _global_id_to_proxy_index_map.end())
                    {
                        _global_id_to_proxy_index_map[f->cell_global_id] = proxy_global_id;
                        _global_proxy_to_global_index_map[ proxy_global_id ] = f->cell_global_id; // save the old global id
                        ++proxy_global_id;
                    }

                    // set out global ID to our proxy id
                    f->cell_global_id = _global_id_to_proxy_index_map[f->cell_global_id];

                    f->owner = _mpirank_to_proxy_mpirank_map[f->owner];
                }

            }

            // ensure the id order is monotonic increasing
            auto perm = sort_permutation(_local_faces, [](mesh_elem const& fa, mesh_elem const& fb)
                                         { return fa->cell_global_id < fb->cell_global_id; });

            _local_faces = apply_permutation_transform(_local_faces, perm);

            for (auto const& name : pars->names)
            {
                _param_data[name] = apply_permutation_transform(_param_data[name], perm);
            }

            auto fname =
                filename_base + ".partition." + (_is_standalone ? "standalone." : "") + std::to_string(mpirank);

            to_hdf5(partition_dir, fname);

            if(output_vtu)
            {
                // the default vtu write param setup isn't going to do what we want so be explicit here
#pragma omp parallel for
                for (int t = 0; t < _local_faces.size(); t++)
                {
                    auto f = _local_faces.at(t);
                    (*f)["global_id"] = f->cell_global_id;
                    (*f)["local_id"] = f->cell_local_id;
                    (*f)["owner"] = f->owner;
                    if (f->is_ghost)
                    {
                        auto& gi = f->get_module_data<ghost_info>("partition_tool");
                        (*f)["ghost_type"] = gi.ghost_type;
                        (*f)["is_ghost"] = f->is_ghost;
                    }
                }

                write_vtu("rank." + std::to_string(mpirank) + ".vtu");
            }

            if(is_subset)
            {
                //undo the proxy changes
                for(auto& f: _local_faces)
                {
                    f->cell_global_id = _global_proxy_to_global_index_map[f->cell_global_id];
                    f->owner = _mpirank_proxy_to_mpirank_map[f->owner];
                }
            }

            pt::ptree m;
            m.put("", (partition_dir / (fname + "_mesh.h5")).string());
            meshes.push_back(std::make_pair("", m));

            pt::ptree p;
            p.put("", (partition_dir / (fname + "_param.h5")).string());
            params.push_back(std::make_pair("", p));
        }
        tree.add_child("meshes", meshes);
        tree.add_child("parameters", params);

        pt::write_json(filename_base + ".np" + std::to_string(_comm_world.size()) + ".partition", tree);
    }

    std::unique_ptr<MeshParameters> read_h5_params(const std::vector<std::string>& param_filenames)
    {
        // Space for extracting the parameter info
        std::unique_ptr<MeshParameters> pars(new MeshParameters);
        for (auto const& param_filename : param_filenames)
        {
            try
            {
                // Turn off the auto-printing when failure occurs so that we can
                // handle the errors appropriately
                Exception::dontPrint();

                // Open an existing file and dataset.
                H5File file(param_filename, H5F_ACC_RDONLY);

                // Extract all of the parameter info (names) from the file
                Group group = file.openGroup("parameters");
                herr_t idx =
                    H5Literate(group.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, group_info, (void*)pars.get());

                for (auto const& name : pars->names)
                {
                    _parameters.insert(name);
                }

                for (auto const& name : pars->names)
                {
                    // hyperslab size that is as large as just the local meshes to read the block
                    hsize_t local_size = static_cast<hsize_t>(_local_faces.size());

                    if (!_is_standalone)
                    {
                        // but we need to ensure we have enough room for local + ghosts
                        _param_data[name].resize(_local_faces.size() + _ghost_faces.size());
                    }
                    else
                    {
                        _param_data[name].resize(_local_faces.size());
                    }
                    hsize_t offset = static_cast<hsize_t>(_local_faces[0]->cell_global_id);

                    DataSet dataset = group.openDataSet(name);
                    DataSpace dataspace = dataset.getSpace();

                    hsize_t nelem;
                    int ndims = dataspace.getSimpleExtentDims(&nelem);

                    dataspace.selectHyperslab(H5S_SELECT_SET, &local_size, &offset);

                    DataSpace memspace(1, &nelem);
                    hsize_t zero = 0;
                    memspace.selectHyperslab(H5S_SELECT_SET, &local_size, &zero);

                    dataset.read(_param_data[name].data(), PredType::NATIVE_DOUBLE, memspace, dataspace);

                    if (!_is_standalone)
                    {
                        // Read parameters for each ghost face individually
                        hsize_t one = 1;
                        // Do NOT do this loop in parallel (internal state of HDF5)
                        for (size_t i = 0; i < _ghost_faces.size(); i++)
                        {
                            auto face = _ghost_faces.at(i);
                            hsize_t global_id = face->cell_global_id;

                            // Position and size in file
                            dataspace.selectHyperslab(H5S_SELECT_SET, &one, &global_id);
                            // Position and size in variable
                            memspace.selectHyperslab(H5S_SELECT_SET, &one, &zero);

                            double value;
                            dataset.read(&value, PredType::NATIVE_DOUBLE, memspace, dataspace);

                            _param_data[name][_local_faces.size() + i] = value;
                        }
                    }
                }
            }
            // catch failure caused by the H5File operations
            catch (FileIException& error)
            {
                error.printErrorStack();
            }
            // catch failure caused by the DataSet operations
            catch (DataSetIException& error)
            {
                error.printErrorStack();
            }


        }
            return pars;
    }

    void to_hdf5(boost::filesystem::path partition_dir, std::string filename_base)
    {

        // update our version
        _version.from_string("3.0.0");

        std::string filename = filename_base + "_mesh.h5";

        //  a global to local that works with the local_faces that has all the ghosts mixed in
        std::map<int, int> global_to_local_faces_index_map;

        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            Exception::dontPrint();

            H5::H5File file((partition_dir / filename).string(), H5F_ACC_TRUNC);
            H5::Group group(file.createGroup("/mesh"));

            // local_faces here is slightly different than anywhere else by, at this  point, it will have ghosts
            // mixed in with it. This is thus all the faces this partition knows about
            hsize_t ntri = _local_faces.size();

            { // Local sizes
                LOG_DEBUG << "Writing Local Sizes.";
                hsize_t npart = _comm_world.size();
                H5::DataSpace dataspace(1, &npart);

                H5::DataSet dataset = file.createDataSet("/mesh/local_sizes", PredType::STD_I32BE, dataspace);
                dataset.write(_local_sizes.data(), PredType::NATIVE_INT);
            }

            { // global elem id
                LOG_DEBUG << "Writting owners";
                H5::DataSpace dataspace(1, &ntri);
                auto owners = std::vector<int>(ntri);

                //                #pragma omp parallel for
                for (size_t i = 0; i < ntri; ++i)
                {
                    owners.at(i) = _local_faces.at(i)->owner;
                }
                H5::DataSet dataset = file.createDataSet("/mesh/owner", PredType::STD_I32BE, dataspace);
                dataset.write(owners.data(), PredType::NATIVE_INT);
            }

            { // global elem id
                LOG_DEBUG << "Writting Global IDs";
                H5::DataSpace dataspace(1, &ntri);
                auto globalIDs = std::vector<int>(ntri);

                //                #pragma omp parallel for
                for (size_t i = 0; i < ntri; ++i)
                {
                    globalIDs.at(i) = _local_faces.at(i)->cell_global_id;
                    global_to_local_faces_index_map[globalIDs.at(i)] = i;
                }
                H5::DataSet dataset = file.createDataSet("/mesh/cell_global_id", PredType::STD_I32BE, dataspace);
                dataset.write(globalIDs.data(), PredType::NATIVE_INT);
            }

            std::set<size_t> vertex_global_id;
            // figure out all of the unique vertexes this subset needs and save their global id
            for (size_t i = 0; i < ntri; ++i)
            {
                auto f = _local_faces.at(i);
                for (size_t j = 0; j < 3; ++j)
                {
                    vertex_global_id.insert(f->vertex(j)->get_id()); // set of the vertexes we will need
                }
            }

            LOG_DEBUG << "This mesh partition has " << vertex_global_id.size() << " unique vertexes";
            // build a list of all the above identified vertexes and set their local id
            std::vector<Delaunay::Vertex_handle> local_vertexes; // vertexes that this part of the mesh needs
            local_vertexes.resize(vertex_global_id.size());
            size_t i = 0;
            for (const auto& itr : vertex_global_id)
            {
                auto vh = vertex(itr);
                vh->set_local_id(i);
                local_vertexes.at(i) = vh;
                ++i;
            }

            { // Which vertices define each face
                H5::DataSpace dataspace(1, &ntri);
                std::vector<std::array<int, 3>> elem(ntri);

#pragma omp parallel for
                for (size_t i = 0; i < ntri; ++i)
                {
                    auto f = _local_faces.at(i);
                    for (size_t j = 0; j < 3; ++j)
                    {
                        elem[i][j] = f->vertex(j)->get_local_id(); // uses a per partition local id
                    }
                }

                H5::DataSet dataset = file.createDataSet("/mesh/elem", elem_t, dataspace);
                dataset.write(elem.data(), elem_t);
            }

            hsize_t nvert = local_vertexes.size();
            { // Vertices
                H5::DataSpace dataspace(1, &nvert);
                std::vector<std::array<double, 3>> vertices(nvert);

#pragma omp parallel for
                for (size_t i = 0; i < nvert; ++i)
                {
                    auto v = local_vertexes.at(i);
                    vertices[i][0] = v->point().x();
                    vertices[i][1] = v->point().y();
                    vertices[i][2] = v->point().z();
                }

                H5::DataSet dataset = file.createDataSet("/mesh/vertex", vertex_t, dataspace);
                dataset.write(vertices.data(), vertex_t);
            }

            { // neighbours
                H5::DataSpace dataspace(1, &ntri);
                std::vector<std::array<int, 3>> neighbor(ntri);

#pragma omp parallel for
                for (size_t i = 0; i < ntri; ++i)
                {
                    auto f = this->face(i);
                    for (size_t j = 0; j < 3; ++j)
                    {
                        auto neigh = f->neighbor(j);
                        if (neigh != nullptr && // we have a neighbour
                            global_to_local_faces_index_map.find(neigh->cell_global_id) !=
                                global_to_local_faces_index_map.end()) // and we've seen this neighbour's global ID
                        {
                            neighbor.at(i)[j] = neigh->cell_global_id;
                        }
                        else
                        {
                            // they may have a triangle neighbour, but we haven't seen it, so  it is outside our ghost
                            // region
                            neighbor.at(i)[j] = -1;
                        }
                    }
                }
                H5::DataSet dataset = file.createDataSet("/mesh/neighbor", neighbor_t, dataspace);
                dataset.write(neighbor.data(), neighbor_t);
            }

            { // ghosts
                H5::DataSpace dataspace(1, &ntri);
                std::vector<int> is_ghost(ntri);
#pragma omp parallel for
                for (size_t ii = 0; ii < ntri; ++ii)
                {
                    // we need to save the following status for ghost faces:
                    // 0 = Not a ghost
                    // 1 = Neighbour ghost
                    // 2 = Radius search ghost
                    auto& gi = _local_faces[ii]->get_module_data<ghost_info>("partition_tool");
                    is_ghost[ii] = gi.ghost_type;
                }
                H5::DataSet dataset = file.createDataSet("/mesh/ghost_type", PredType::STD_I32BE, dataspace);
                dataset.write(is_ghost.data(), PredType::NATIVE_INT);
            }

            // Ensure the proj4 string can fit in the HDF5 data type
            if (_srs_wkt.length() >= 256)
            {
                BOOST_THROW_EXCEPTION(config_error() << errstr_info("Proj4 string needs to be < 256. Length: " +
                                                                    std::to_string(_srs_wkt.length())));
            }
            { // Write the proj4
                H5::DataSpace dataspace(1, &proj4_dims);
                H5::Attribute attribute = file.createAttribute("/mesh/proj4", proj4_t, dataspace);
                attribute.write(proj4_t, _srs_wkt);
            }

            { // Write the partition method
                H5::DataSpace dataspace(1, &partition_type_dims);
                H5::Attribute attribute = file.createAttribute("/mesh/partition_method", partition_type_t, dataspace);
                attribute.write(partition_type_t, _partition_method);
            }

            { // Write the is_geographic
                H5::DataSpace dataspace(1, &geographic_dims);
                H5::Attribute attribute =
                    file.createAttribute("/mesh/is_geographic", PredType::NATIVE_HBOOL, dataspace);
                attribute.write(PredType::NATIVE_HBOOL, &_is_geographic);
            }

            { // Write that this is part of a partitioned mesh
                _is_partition = true;
                if (_is_standalone)
                    _is_partition = false;

                H5::DataSpace dataspace(1, &partition_dims);
                H5::Attribute attribute = file.createAttribute("/mesh/is_partition", PredType::NATIVE_HBOOL, dataspace);
                attribute.write(PredType::NATIVE_HBOOL, &_is_partition);
            }

            { // Write the mesh version
                H5::DataSpace dataspace(1, &meshver_dims);
                H5::Attribute attribute = file.createAttribute("/mesh/version", meshver_t, dataspace);
                attribute.write(meshver_t, _version.to_string());
            }

        } // end of try block

        // catch failure caused by the H5File operations
        catch (const FileIException& error)
        {
            error.printErrorStack();
        }

        // catch failure caused by the DataSet operations
        catch (const DataSetIException& error)
        {
            error.printErrorStack();
        }

        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////
        // Parameters
        ////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////

        std::string par_filename = filename_base + "_param.h5";

        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            Exception::dontPrint();

            H5::H5File file((partition_dir / par_filename).string(), H5F_ACC_TRUNC);
            H5::Group group(file.createGroup("/parameters"));

            hsize_t ntri = _local_faces.size();

            for (auto& par_iter : _param_data)
            {
                const auto& p = par_iter.first;

                H5::DataSpace dataspace(1, &ntri);
                std::string par_location = "/parameters/" + p;
                H5::DataSet dataset = file.createDataSet(par_location, PredType::NATIVE_DOUBLE, dataspace);

                dataset.write(_param_data[p].data(), PredType::NATIVE_DOUBLE);
            }

        } // end try block
        // catch failure caused by the H5File operations
        catch (FileIException& error)
        {
            error.printErrorStack();
        }
        // catch failure caused by the DataSet operations
        catch (DataSetIException& error)
        {
            error.printErrorStack();
        }
    }

    class comm_world
    {
      public:
        int rank() { return _rank; }

        int size() { return _size; }

        int _size;
        int _rank;
    };

    comm_world _comm_world;

    std::map<std::string, std::vector<double>> _param_data;

    bool _is_partition; // the flag to write to the hdf5 file

    bool _is_standalone; // are we in standalone mode

    struct ghost_info : public face_info
    {
        // Note the type of ghost we are:
        // None = not a ghost
        // Neigh = Neighbour ghost
        // Dist = Distance search ghost
        enum GHOST_TYPE
        {
            NONE,
            NEIGH,
            DIST
        };
        GHOST_TYPE ghost_type = GHOST_TYPE::NONE;
    };
};

int main(int argc, char* argv[])
{
    BOOST_LOG_FUNCTION()

    std::string mesh_filename;
    size_t MPI_ranks = 2;
    double max_ghost_distance = 100;

    int standalone_rank = -1;
    bool output_vtu=false;

    po::options_description desc("Allowed options.");
    desc.add_options()("help", "This message")
        ("mesh-file,m", po::value<std::string>(&mesh_filename), "Mesh file")(
            "param-file,p", po::value<std::vector<std::string>>(),"Parameter file")(
            "max-ghost-distance,g", po::value<double>(&max_ghost_distance), "Maximum ghost distance")(
            "valid-ranks,v", po::value<std::vector<int>>(),"Only output these ranks")(
        "standalone,s", po::value<int>(&standalone_rank),
        "Write the paritioned mesh so-as to be loadable standalone. Advanced debugging feature. Arg is rank to do")(
            "write-vtu", po::bool_switch(&output_vtu),"Output the rank partitions as vtu")(
        "mpi-ranks", po::value<size_t>(&MPI_ranks), "Number of MPI ranks to partition for");

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
    po::notify(vm);

    auto is_json_mesh = boost::filesystem::path(mesh_filename).extension().string() == ".mesh";

    if (is_json_mesh)
    {
        LOG_WARNING << "The input mesh is in json format. It MUST be converted to hdf5 before it can be partitioned.\n"
                       " Once the hdf5 conversion is done, this tool will re-run with the h5 mesh as input to fully "
                       "partition the mesh.";
    }

    if (!vm.count("mesh-file"))
    {
        LOG_ERROR << "Mesh file required";
        exit(-1);
    }
    if (!vm.count("param-file"))
    {
        LOG_ERROR << "Param file required";
        exit(-1);
    }
    if (!vm.count("mpi-ranks") && !is_json_mesh)
    {
        LOG_ERROR << "MPI ranks required";
        exit(-1);
    }

    if (vm.count("mpi-ranks") && is_json_mesh)
    {
        LOG_WARNING << "MPI ranks will be ignored for the json -> h5 conversion. ";
    }

    if (!vm.count("max-ghost-distance"))
    {
        LOG_WARNING << "Using default max ghost distance of 100 m";
    }

    std::vector<int> valid_ranks;
    if(vm.count("valid-ranks"))
    {
        for (auto& itr : vm["valid-ranks"].as<std::vector<int>>())
        {
            valid_ranks.push_back(itr);
        }
    }

    if(valid_ranks.empty())
    {
        valid_ranks.resize(MPI_ranks);
        std::iota(valid_ranks.begin(), valid_ranks.end(), 0);
    }


    std::vector<std::string> param_filenames;
    for (auto& itr : vm["param-file"].as<std::vector<std::string>>())
    {
        param_filenames.push_back(itr);
    }

    if (vm.count("standalone"))
    {
        LOG_WARNING << "Standalone option enabled. This will not write ghost faces and is intended to write a specific "
                       "rank as a standalone mesh for debugging.";
    }

    if (MPI_ranks <= 1)
    {
        LOG_ERROR << "Requires ranks >1";
        exit(-1);
    }

    try
    {

        {
            preprocessingTriangulation* tri = new preprocessingTriangulation();
            tri->from_file_and_partition(mesh_filename, param_filenames, MPI_ranks, max_ghost_distance,
                                         standalone_rank, valid_ranks, output_vtu);
        }

        // We have input as json and asked for mpi-ranks, so rerun the tool
        if (is_json_mesh && vm.count("mpi-ranks"))
        {
            LOG_DEBUG << "Partitioning h5 mesh";

            auto h5_mesh_name = boost::filesystem::path(mesh_filename).stem().string() + "_mesh.h5";
            std::vector<std::string> h5_param_name = {boost::filesystem::path(mesh_filename).stem().string() +
                                                      "_param.h5"};

            preprocessingTriangulation* tri = new preprocessingTriangulation();
            tri->from_file_and_partition(h5_mesh_name, h5_param_name, MPI_ranks, max_ghost_distance, standalone_rank, valid_ranks, output_vtu);

        }
    }
    catch (const exception_base& e)
    {
        LOG_ERROR << boost::diagnostic_information(e);
    }

    LOG_DEBUG << "Done";
}
