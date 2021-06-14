#include "H5Cpp.h"
#include "logger.hpp"
#include "triangulation.hpp"
using namespace H5;

class preprocessingTriangulation : public triangulation
{
  public:
    void from_hdf5_and_partition(const std::string& mesh_filename, int MPI_ranks)
    {
        std::vector<int> permutation;
        std::vector<Point_2> center_points;

        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            Exception::dontPrint();

            // Open an existing file and dataset.
            H5File file(mesh_filename, H5F_ACC_RDONLY);

            std::vector<std::array<double, 3>> vertex;
            std::vector<std::array<int, 3>> elem;
            std::vector<std::array<int, 3>> neigh;

            {
                H5::DataSet dataset = file.openDataSet("/mesh/cell_global_id");
                H5::DataSpace dataspace = dataset.getSpace();
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

                    Vertex_handle Vh = this->create_vertex();
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

                    auto face = this->create_face(vert1, vert2, vert3);
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
                LOG_DEBUG << "Created a mesh with " << this->size_faces() << " triangles";
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

        //        partition_mesh();
        LOG_DEBUG << "Partitioning mesh";
        _num_global_faces = _faces.size();


        std::cout << MPI_ranks << std::endl;
        for (int mpirank = 0; mpirank < MPI_ranks; mpirank++)
        {
            int my_rank = mpirank;

            // Set up so that all processors know how 'big' all other processors are

            _num_faces_in_partition.resize(MPI_ranks, _num_global_faces / MPI_ranks);
            for (unsigned int i = 0; i < _num_global_faces % MPI_ranks; ++i)
            {
                _num_faces_in_partition[i]++;
                LOG_DEBUG <<_num_faces_in_partition[i] ;
            }

            // each processor only knows its own start and end indices
            size_t face_start_idx = 0;
            size_t face_end_idx = _num_faces_in_partition[0] - 1;
            for (int i = 1; i <= mpirank; ++i)
            {
                face_start_idx += _num_faces_in_partition[i - 1];
                face_end_idx += _num_faces_in_partition[i];
            }
            // Store in the public members
            global_cell_start_idx = face_start_idx;
            global_cell_end_idx = face_end_idx;

            // Set size of vector containing locally owned faces
            _local_faces.resize(_num_faces_in_partition[mpirank]);

            _global_IDs.resize(_local_faces.size());
            // Loop can't be parallel due to modifying map
            for (size_t local_ind = 0; local_ind < _local_faces.size(); ++local_ind)
            {
                _global_IDs[local_ind] = face_start_idx + local_ind;
                _global_to_locally_owned_index_map[_global_IDs[local_ind]] = local_ind;

                _faces.at(_global_IDs[local_ind])->is_ghost = false;
                _faces.at(_global_IDs[local_ind])->owner = my_rank;
                _faces.at(_global_IDs[local_ind])->cell_local_id = local_ind;
                _local_faces[local_ind] = _faces.at(_global_IDs[local_ind]);
            }

            LOG_DEBUG << "MPI Process " << my_rank << ": start " << face_start_idx << ", end " << face_end_idx
                      << ", number " << _local_faces.size();
        }

//        _num_faces = _local_faces.size();
//        determine_local_boundary_faces();
//        determine_process_ghost_faces_nearest_neighbors();
//
//        // Region
//        // TODO: Need to auto-determine how far to look based on module setups
//        determine_process_ghost_faces_by_distance(100.0);


    }
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
