#include "H5Cpp.h"
#include "logger.hpp"
#include "triangulation.hpp"
using namespace H5;

class preprocessingTriangulation : public triangulation
{
  public:
    void partition(int mpi_rank)
    { // Set up so that all processors know how 'big' all other processors are
        _num_faces_in_partition.resize(_comm_world.size(), _num_global_faces / _comm_world.size());
        for (unsigned int i = 0; i < _num_global_faces % _comm_world.size(); ++i)
        {
            _num_faces_in_partition[i]++;
        }

            // each processor only knows its own start and end indices
            size_t face_start_idx = 0;
            size_t face_end_idx = _num_faces_in_partition[0] - 1;
            for (int i = 1; i <= _comm_world.rank(); ++i)
            {
                face_start_idx += _num_faces_in_partition[i - 1];
                face_end_idx += _num_faces_in_partition[i];
            }
            // Store in the public members
            global_cell_start_idx = face_start_idx;
            global_cell_end_idx = face_end_idx;

            // Set size of vector containing locally owned faces
            _local_faces.resize(_num_faces_in_partition[_comm_world.rank()]);

            _global_IDs.resize(_local_faces.size());
            // Loop can't be parallel due to modifying map
            for (size_t local_ind = 0; local_ind < _local_faces.size(); ++local_ind)
            {
                _global_IDs[local_ind] = face_start_idx + local_ind;
                _global_to_locally_owned_index_map[_global_IDs[local_ind]] = local_ind;

                _faces.at(_global_IDs[local_ind])->is_ghost = false;
                _faces.at(_global_IDs[local_ind])->owner = mpi_rank;
                _faces.at(_global_IDs[local_ind])->cell_local_id = local_ind;
                _local_faces[local_ind] = _faces.at(_global_IDs[local_ind]);
            }

            LOG_DEBUG << "MPI Process " << mpi_rank << ": start " << face_start_idx << ", end " << face_end_idx
                      << ", number " << _local_faces.size();
    }
    void read_h5(const std::string& mesh_filename)
    {
        std::vector<int> permutation;
        std::__1::vector<Point_2> center_points;

        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            Exception::dontPrint();

            // Open an existing file and dataset.
            H5File file(mesh_filename, H5F_ACC_RDONLY);

            std::__1::vector<std::array<double, 3>> vertex;
            std::__1::vector<std::array<int, 3>> elem;
            std::__1::vector<std::array<int, 3>> neigh;

            {
                DataSet dataset = file.openDataSet("/mesh/cell_global_id");
                DataSpace dataspace = dataset.getSpace();
                hsize_t nelem;
                int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);
                permutation.resize(nelem);
                dataset.read(permutation.data(), PredType::NATIVE_INT);
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
                    //                    face->_domain = this;

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

                    // LOG_DEBUG << neigh[i][0] << " " << neigh[i][1] << " " << neigh[i][2];

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
    void determine_local_boundary_faces()
    {
        /*
          - Store handles to boundary faces on the locally owned process
          - Also store boolean value "is_global_boundary"
        */

        using th_safe_multicontainer_type = std::vector< std::pair<mesh_elem,bool> >[];


        // Need to ensure we're starting from nothing?
        assert( _boundary_faces.size() == 0 );

        LOG_DEBUG << "Determining local boundary faces";

        std::unique_ptr< th_safe_multicontainer_type > th_local_boundary_faces;

#pragma omp parallel
        {
            // We want an array of vectors, so that OMP threads can increment them
            // separately, then join them afterwards
#pragma omp single
            {
                th_local_boundary_faces = std::make_unique< th_safe_multicontainer_type >(omp_get_num_threads());
            }
#pragma omp for
            for(size_t face_index=0; face_index< _local_faces.size(); ++face_index)
            {
                // face_index is a local index... get the face handle
                auto face = _local_faces.at(face_index);

                int num_owned_neighbors = 0;
                for (int neigh_index = 0; neigh_index < 3; ++neigh_index)
                {

                    auto neigh = face->neighbor(neigh_index);

                    // Test status of neighbor
                    if (neigh == nullptr)
                    {
                        th_local_boundary_faces[omp_get_thread_num()].push_back(std::make_pair(face,true));
                        num_owned_neighbors=3; // set this to avoid triggering the post-loop if statement
                        break;
                    } else
                    {
                        if (neigh->is_ghost == false)
                        {
                            num_owned_neighbors++;
                        }
                    }
                }

                // If we don't own 3 neighbors, we are a local, but not a global boundary face
                if( num_owned_neighbors<3 ) {
                    th_local_boundary_faces[omp_get_thread_num()].push_back(std::make_pair(face,false));
                }
            }

            // Join the vectors via a single thread in t operations
            //  NOTE future optimizations:
            //   - reserve space for insertions into _boundary_faces
            //   - can be done recursively in log2(t) operations
#pragma omp single
            {
                for(int thread_idx=0;thread_idx<omp_get_num_threads();++thread_idx) {
                    _boundary_faces.insert(_boundary_faces.end(),
                                           th_local_boundary_faces[thread_idx].begin(), th_local_boundary_faces[thread_idx].end());
                }
            }
        }

        // Some log debug output to see how many boundary faces on each
        LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _boundary_faces.size() << " boundary faces.";

        LOG_DEBUG << "MPI Process " << _comm_world.rank() << " _faces.size(): " << _faces.size() << " _local_faces.size(): " << _local_faces.size();


    }
    void determine_process_ghost_faces_nearest_neighbors()
    {
        // NOTE that this algorithm is not implemented for multithread
        // - multithread can be implemented similarly to determine_local_boundary_faces
        //   - not a priority since number of boundary faces should be small

        // Ensure that the local boundary faces have been determined, but ghost nearest neighbors have not been set
        assert( _boundary_faces.size() != 0 );
        assert( _ghost_neighbors.size() == 0 );

        LOG_DEBUG << "Determining ghost region info";

        // Vector for append speed
        std::vector< mesh_elem > ghosted_boundary_nearest_neighbors;

        for(size_t face_index=0; face_index< _boundary_faces.size(); ++face_index)
        {
            // face_index is a local index... get the face handle
            auto face = _boundary_faces.at(face_index).first;
            // append the ghosted nearest neighbors
            for(int i = 0; i < 3; ++i)
            {
                auto neigh = face->neighbor(i);
                if(neigh != nullptr && neigh->is_ghost)
                    ghosted_boundary_nearest_neighbors.push_back(neigh);
            }
        }

        // Convert to a set to remove duplicates
        std::unordered_set<mesh_elem> tmp_set(std::begin(ghosted_boundary_nearest_neighbors),
                                              std::end(ghosted_boundary_nearest_neighbors));
        // Convert the set to a vector
        _ghost_neighbors.insert(std::end(_ghost_neighbors),
                                std::begin(tmp_set),std::end(tmp_set));
        // Sort the ghost neighbors
        // NOTE:
        // - sorting this vector by cell_global_id effectively partitions the ghost neighbors to be contiguous in communication partners
        std::sort(_ghost_neighbors.begin(),_ghost_neighbors.end(),
                  [&](const auto& a, const auto& b)
                  {
                    return a->cell_global_id < b->cell_global_id;
                  });


        LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _ghost_neighbors.size() << " ghosted nearest neighbors.";


        // Determine the owners of the ghost faces (for communication setup)
        _ghost_neighbor_owners.resize(_ghost_neighbors.size());
        int start_index=0;
        int prev_owner;
        int num_partners=0;
        // Construct ghost region ownership info
        for(size_t i=0; i<_ghost_neighbors.size(); ++i)
        {
            // index type needs to match type of elements of _num_faces_in_partition
            int global_ind = static_cast<int>(_ghost_neighbors[i]->cell_global_id);
            _ghost_neighbor_owners[i] = determine_owner_of_global_index(global_ind,
                                                                        _num_faces_in_partition);
            // on first it, no value of prev_owner exists
            if(i==0) prev_owner = _ghost_neighbor_owners[i];
            // if owner different from last iteration, store prev segment's ownership info
            if (prev_owner != _ghost_neighbor_owners[i])
            {
                num_partners++;
                _comm_partner_ownership[prev_owner] = std::make_pair(start_index, i-start_index);
                start_index=i;
            }

            if (i ==_ghost_neighbors.size()-1) {
                _comm_partner_ownership[prev_owner] = std::make_pair(start_index, i-start_index+1);
            }

            // prep prev_owner for next iteration
            prev_owner=_ghost_neighbor_owners[i];
        }

        for(size_t i=0; i<_ghost_neighbors.size(); ++i)
        {
            _global_index_to_local_ghost_map[_ghost_neighbors[i]->cell_global_id] = static_cast<int>(i);
        }


        LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _comm_partner_ownership.size() << " communication partners:";
        for (auto it : _comm_partner_ownership)
        {
            LOG_DEBUG << "  # Process " << _comm_world.rank() << " partner " << it.first << " owns local ghost indices " << it.second.first << " to " << it.second.first+it.second.second-1;
        }


    }

    void from_hdf5_and_partition(const std::string& mesh_filename, int MPI_ranks)
    {
        read_h5(mesh_filename);

        LOG_DEBUG << "Partitioning mesh";
        _num_global_faces = _faces.size();

        _comm_world._size = 2;

        std::cout << MPI_ranks << std::endl;
        for (int mpirank = 0; mpirank < MPI_ranks; mpirank++)
        {
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
            _global_IDs = {};

#pragma omp parallel for
            for(size_t i = 0; i < _faces.size(); ++i)
            {
                _faces.at(i)->is_ghost = true; // default state, switches to false later
            }

            partition(_comm_world._rank);

            _num_faces = _local_faces.size();

            determine_local_boundary_faces();
            determine_process_ghost_faces_nearest_neighbors();
        }



//        // Region
//        // TODO: Need to auto-determine how far to look based on module setups
//        determine_process_ghost_faces_by_distance(100.0);


    }

    class comm_world
    {
      public:

        int rank()
        {
            return _rank;
        }

        int size()
        {
            return _size;
        }

        int _size;
        int _rank;
    };

    comm_world _comm_world;
};

int main()
{
    BOOST_LOG_FUNCTION();

    std::string mesh_filename;

    int MPI_ranks = 2;

    mesh_filename = "/Users/chris/Documents/science/model_runs/benchmark_problems/granger_pbsm_synthetic/test_mesh.h5";

    preprocessingTriangulation* tri = new preprocessingTriangulation();
    tri->from_hdf5_and_partition(mesh_filename,MPI_ranks);
}
