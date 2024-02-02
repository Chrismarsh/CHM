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



#include "triangulation.hpp"

triangulation::triangulation()
{

    _num_faces = 0;
    _vtk_unstructuredGrid = nullptr;
    _is_geographic = false;
    _UTM_zone = 0;
    _terrain_deformed=false;
    _min_z =  999999;
    _max_z = -999999;

#ifdef USE_SPARSEHASH
    data.set_empty_key("");
    vectors.set_empty_key("");
    vertex_data.set_empty_key("");

#endif

/*
  Datatypes for reading/writing HDF5 files
*/
    // Array datatype for vertices
    vertex_dims=3;
    vertex_t = H5::ArrayType(PredType::NATIVE_DOUBLE,1,&vertex_dims);

    // Array datatype for vertices defining faces
    elem_dims=3;
    elem_t = H5::ArrayType(PredType::NATIVE_INT,1,&elem_dims);

    // Array datatype for neighbors
    neighbor_dims=3;
    neighbor_t = H5::ArrayType(PredType::NATIVE_INT,1,&neighbor_dims);

    // Array datatype for proj4 string
    proj4_dims=1;
    proj4_t = H5::StrType(PredType::C_S1,256);

    meshver_dims=1;
    meshver_t = H5::StrType(PredType::C_S1,256);

    // Array datatype for partition method string
    partition_type_dims=1;
    partition_type_t = H5::StrType(PredType::C_S1,256);
    _partition_method = "";

    // Array datatype for is_geographic
    geographic_dims=1;
    partition_dims=1;

}



triangulation::~triangulation()
{

}

std::string triangulation::proj4()
{
    return _srs_wkt;
}

void triangulation::write_param_to_vtu(bool write_param)
{
    _write_parameters_to_vtu = write_param;
}

void triangulation::write_ghost_neighbors_to_vtu(bool write_ghost_neighbors)
{
    _write_ghost_neighbors_to_vtu = write_ghost_neighbors;
}

std::set<std::string> triangulation::parameters()
{
    return _parameters;
}

bool triangulation::is_geographic()
{
    return _is_geographic;
}
size_t triangulation::size_faces()
{
    return _num_faces;
}
size_t triangulation::size_global_faces()
{
    return _num_global_faces;
}

size_t triangulation::size_vertex()
{
    return _num_vertex;
}

mesh_elem triangulation::locate_face(double x, double y)
{
    Point_2 query(x,y);
    return locate_face(query);
}

mesh_elem triangulation::locate_face(Point_2 query)
{
    // In some limited cases triangles can be arranged in a circle around a central vertex
    // In this case, there could be 360/21.5 = 16 and change triangles.
    // So just grab the nearest 17 triangles, one of these will hold the point we need
    K_neighbor_search search(*dD_tree, query, 17);

    //check if the closest is what we wanted
    for(auto itr: search)
    {
      auto f = boost::get<1>(itr.first); //grab the triangle from the iterator
      if(!f->is_ghost &&
          f->contains(query.x(),query.y()))
        return boost::get<1>(itr.first);
    }


    //we tried....
    return nullptr;
}


std::vector<mesh_elem > triangulation::find_faces_in_radius(Point_2 center, double radius) const
{
    // define exact circular range query  (fuzziness=0)
    Fuzzy_circle exact_range(center, radius);

    std::vector<boost::tuple<K::Point_2, mesh_elem > > result;
    dD_tree->search(std::back_inserter(result), exact_range);

    std::vector< mesh_elem > faces;

    for (auto& itr : result)
    {
        faces.push_back( boost::get<1>(itr));
    }
    return faces;
}


std::vector< mesh_elem > triangulation::find_faces_in_radius(double x, double y, double radius) const
{
    Point_2 query(x,y);
    return find_faces_in_radius(query, radius);
}



mesh_elem triangulation::find_closest_face(Point_2 query) const
{
    K_neighbor_search search(*dD_tree, query, 1);

    //return the first face only
    auto it = search.begin();
    return  boost::get<1>(it->first);

}
mesh_elem triangulation::find_closest_face(double x, double y) const
{

    Point_2 query(x,y);

    return find_closest_face(query);

}


void triangulation::serialize_parameter(std::string output_path, std::string parameter)
{

    pt::ptree out;
    pt::ptree subtree;

    for(auto& face : _faces)
    {
        pt::ptree item;
        item.put_value(face->parameter(parameter));
        subtree.push_back(std::make_pair("",item));
    }
    out.put_child( pt::ptree::path_type(parameter),subtree);

    pt::write_json(output_path,out);

}
void triangulation::from_json(pt::ptree &mesh)
{

    _version.from_string(mesh.get<std::string>("mesh.version",""));

    if(!_version.mesh_ver_meets_min_json())
    {
        CHM_THROW_EXCEPTION(mesh_error, "JSON mesh version to too old");
    }

    size_t nvertex_toread = mesh.get<size_t>("mesh.nvertex");
    SPDLOG_DEBUG("Reading in #vertex={}",nvertex_toread);
    size_t i=0;

    //paraview struggles with lat/long as it doesn't seem to have the accuracy. So we need to scale up lat-long.
    int is_geographic = mesh.get<int>("mesh.is_geographic");
    if( is_geographic == 1)
    {
        _is_geographic = true;
    }
    _srs_wkt = mesh.get<std::string>("mesh.proj4","");

    _partition_method = mesh.get<std::string>("mesh.partition_method","");

    if(_srs_wkt == "")
    {
        CHM_THROW_EXCEPTION(config_error, "proj4 field in .mesh file is empty!");
    }

    for (auto &itr : mesh.get_child("mesh.vertex"))
    {
        std::vector<double> vertex;
        //iterate over the vertex triples
        for(auto& jtr: itr.second)
        {
            vertex.push_back(jtr.second.get_value<double>());
        }
        Point_3 pt( vertex[0], vertex[1], vertex[2]);

        _max_z = std::max(_max_z,vertex[2]);
        _min_z = std::min(_min_z,vertex[2]);

        _bounding_box.x_max = std::max(_bounding_box.x_max, vertex[0]);
        _bounding_box.x_min = std::min(_bounding_box.x_min, vertex[0]);

        _bounding_box.y_max = std::max(_bounding_box.y_max, vertex[1]);
        _bounding_box.y_min = std::min(_bounding_box.y_min, vertex[1]);


        Vertex_handle Vh = this->create_vertex();
        Vh->set_point(pt);
        Vh->set_id(i);
        _vertexes.push_back(Vh);
        i++;
    }
    _num_vertex = this->number_of_vertices();
    SPDLOG_DEBUG("# nodes created = {}",_num_vertex);

    if( this->number_of_vertices() != nvertex_toread)
    {
        CHM_THROW_EXCEPTION(config_error,
                    "Expected: " + std::to_string(nvertex_toread) + " vertex, got: " + std::to_string(this->number_of_vertices()));
    }

    //read in faces
    size_t num_elem = mesh.get<int>("mesh.nelem");
    SPDLOG_DEBUG("Reading in #elem = {}",num_elem);
    //set our mesh dimensions
    this->set_dimension(2);


    //vectors to hold the center of a face to generate the spatial search tree
    std::vector<Point_2> center_points;

    i = 0;
    size_t cid = 0;
    for (auto &itr : mesh.get_child("mesh.elem"))
    {
        std::vector<int> items;
        //iterate over the vertex triples
        for(auto& jtr: itr.second)
        {
            items.push_back(jtr.second.get_value<size_t>());
        }
        auto vert1 = _vertexes.at(items[0]); //0 indexing
        auto vert2 = _vertexes.at(items[1]);
        auto vert3 = _vertexes.at(items[2]);

        auto face = this->create_face(vert1,vert2,vert3);
        face->cell_global_id = cid++;
        face->cell_local_id = face->cell_global_id;

        _global_IDs.push_back(face->cell_global_id);

        if( is_geographic == 1)
        {
            face->_is_geographic = true;
        }
        face->_debug_ID= --i; //all ids will be negative starting at -1. Named ids (for output) will be positive starting at 0
        face->_debug_name= std::to_string(i);
        face->_domain = this;

        vert1->set_face(face);
        vert2->set_face(face);
        vert3->set_face(face);


        _faces.push_back(face);
    }



    _num_faces = this->number_of_faces();

    SPDLOG_DEBUG("Created a mesh with {} triangles", this->size_faces());

    if( this->number_of_faces() != num_elem)
    {
        CHM_THROW_EXCEPTION(config_error,
                "Expected: " + std::to_string(num_elem) + " elems, got: " + std::to_string(this->size_faces()));
    }

    SPDLOG_DEBUG("Building face neighbors");
    i=0;
    int nelem = mesh.get<int>("mesh.nelem"); // what we are expecting to see, 0 indexed
    for (auto &itr : mesh.get_child("mesh.neigh"))
    {
        std::vector<int> items;
        //iterate over the vertex triples
        for(auto& jtr: itr.second)
        {
            items.push_back(jtr.second.get_value<int>());
        }

        auto face = _faces.at(i);

        if(    items[0] > nelem
            || items[1] > nelem
            || items[2] > nelem)
        {
            CHM_THROW_EXCEPTION(config_error, "Face " + std::to_string(i) + " has out of bound neighbors.");
        }
        //-1 is now the no neighbor value
        Face_handle face0 =  items[0] != -1 ?_faces.at( items[0] ) : nullptr; //Face_handle()
        Face_handle face1 =  items[1] != -1 ?_faces.at( items[1] ) : nullptr;
        Face_handle face2 =  items[2] != -1 ?_faces.at( items[2] ) : nullptr;

        face->set_neighbors(face0,face1,face2);

        i++;
    }
    //loop over parameter name


    try
    {
        // build up the entire list of parameters so we can use this to init the per-face parameter
        // storage later

        std::set<std::string> blacklist; // holds any parameters that have 0 length, eg "area": [],
        for (auto &itr : mesh.get_child("parameters"))
        {
            // we could have an item like this
            // "area": [],
            // and we need to ensure we *don't* load those
            std::string name = itr.first.data();
            size_t i=0;
            for (auto &jtr : itr.second)
            {
                //just count the first couple items, make sure it's non zero
                if(i > 1)
                    break;
                ++i;
            }

            if(i == 0)
            {
                blacklist.insert(name);
                SPDLOG_WARN("Parameter " + name + " is zero length and will be ignored.");
            } else {
                _parameters.insert(name);
            }

        }

        // init the storage
#pragma omp parallel for
        for (size_t i = 0; i < size_faces(); i++)
        {
             _faces.at(i)->init_parameters(_parameters);
        }


        for (auto &itr : mesh.get_child("parameters"))
        {
            i = 0; // reset evertime we get a new parameter set
            auto name = itr.first.data();

            if(blacklist.find(name) != blacklist.end())
                continue; // skip blacklisted ones, as we don't want to use these

            SPDLOG_DEBUG("Applying parameter: {}",name);

            for (auto &jtr : itr.second)
            {
                try
                {
                    auto face = _faces.at(i);
                    auto value = jtr.second.get_value<double>();
                    face->parameter(name) = value;
                    i++;
                }catch(std::out_of_range& e)
                {
                    SPDLOG_ERROR("Something is wrong with the parameter file. There are more parameter elements than triangulation elements" );
                    CHM_THROW_EXCEPTION(mesh_error, "Something is wrong with the parameter file. There are more parameter elements than triangulation elements");
                }

            }
        }
    }catch(pt::ptree_bad_path& e)
    {
        // we don't have this section, no worries
        // but we still need to build up the face storage as we may have parameters from a module
        // init the storage, which builds the mphf
#pragma omp parallel for
        for (size_t i = 0; i < size_faces(); i++)
        {
             _faces.at(i)->init_parameters(_parameters);
        }

    }

    std::set<std::string> ics;
    try
    {
        for (auto &itr : mesh.get_child("initial_conditions"))
        {
            i=0; // reset evertime we get a new parameter set
            auto name = itr.first.data();
            SPDLOG_DEBUG("Applying IC: {}",name);
            for (auto &jtr : itr.second)
            {
                auto face = _faces.at(i);
                double value = jtr.second.get_value<double>();
                face->set_initial_condition(name,value);
                i++;
            }

            ics.insert(name);
        }
    }catch(pt::ptree_bad_path& e)
    {
        // we don't have this section, no worries
    }
    _num_faces = this->number_of_faces();
    _num_vertex = this->number_of_vertices();

    // Get local sizes for each rank
    // If not available, use the old "balanced" setting
    try
    {
      for (auto itr : mesh.get_child("mesh.local_size"))
      {
	    _local_sizes.push_back(itr.second.get_value<size_t>());
      }
    }catch(pt::ptree_bad_path& e)
    {
        // If not, just ignore the exception
        SPDLOG_DEBUG("No face local sizes.");
    }

    // Permute the faces if they have explicit IDs set in the mesh file
    try
    {
      std::vector<size_t> permutation;
      for (auto itr : mesh.get_child("mesh.cell_global_id"))
      {
	    permutation.push_back(itr.second.get_value<size_t>());
      }
      reorder_faces(permutation);

    }catch(pt::ptree_bad_path& e)
    {
        // If not, just ignore the exception
        SPDLOG_DEBUG("No face permutation.");
    }

    // Don't actually want to partition the mesh. Use the non-MPI one
    partition_mesh_nonMPI(_faces.size());

    _build_dDtree();


#pragma omp parallel for
    for (size_t i = 0; i < _faces.size(); i++)
    {
        auto f = _faces.at(i); // ensure we access like this as json mode is generally non MPI unless we are partitioning
        //init these
        f->slope();
        f->aspect();
        f->center();
        f->normal();
    }


}

void triangulation::to_hdf5(std::string filename_base)
{

  std::string filename = filename_base + "_mesh.h5";

  try {
    // Turn off the auto-printing when failure occurs so that we can
    // handle the errors appropriately
    Exception::dontPrint();

    H5::H5File file(filename,H5F_ACC_TRUNC);
    H5::Group group(file.createGroup("/mesh"));

    hsize_t ntri= size_global_faces();
    hsize_t nvert= size_vertex();

    { // Local sizes
      SPDLOG_DEBUG("Writing Local Sizes.");
      hsize_t npart = _local_sizes.size();
      H5::DataSpace dataspace(1, &npart);

      H5::DataSet dataset = file.createDataSet("/mesh/local_sizes", PredType::STD_I32BE, dataspace);
      dataset.write(_local_sizes.data(), PredType::NATIVE_INT);
    }

    {  // Faces
      H5::DataSpace dataspace(1, &ntri);
      auto globalIDs = get_global_IDs();
      H5::DataSet dataset = file.createDataSet("/mesh/cell_global_id", PredType::STD_I32BE, dataspace);
      dataset.write(globalIDs.data(), PredType::NATIVE_INT);
    }

    { // Vertices
      H5::DataSpace dataspace(1, &nvert);
      std::vector<std::array<double,3>> vertices(nvert);

#pragma omp parallel for
      for (size_t i = 0; i < nvert; ++i)
	{
	  auto v = vertex(i);
	  vertices[i][0] = v->point().x();
	  vertices[i][1] = v->point().y();
	  vertices[i][2] = v->point().z();
	}

      H5::DataSet dataset = file.createDataSet("/mesh/vertex", vertex_t, dataspace);
      dataset.write(vertices.data(), vertex_t);
    }

    {  // Which vertices define each face
      H5::DataSpace dataspace(1, &ntri);
      std::vector<std::array<int,3>> elem(ntri);

#pragma omp parallel for
      for (size_t i = 0; i < ntri; ++i)
	{
	  auto f = this->face(i);
	  for (size_t j = 0; j < 3; ++j) {
	    elem[i][j] = f->vertex(j)->get_id();
	  }
	}

      H5::DataSet dataset = file.createDataSet("/mesh/elem", elem_t, dataspace);
      dataset.write(elem.data(), elem_t);
    }

    {  // Which vertices define each face
      H5::DataSpace dataspace(1, &ntri);
      std::vector<std::array<int,3>> neighbor(ntri);

#pragma omp parallel for
      for (size_t i = 0; i < ntri; ++i)
	{
	  auto f = this->face(i);
	  for (size_t j = 0; j < 3; ++j) {
	    auto neigh = f->neighbor(j);
	    if(neigh != nullptr) {
	      neighbor[i][j] = neigh->cell_global_id;
	    } else {
	      neighbor[i][j] = -1;
	    }
	  }
	}
      H5::DataSet dataset = file.createDataSet("/mesh/neighbor", neighbor_t, dataspace);
      dataset.write(neighbor.data(), neighbor_t);
    }

    // Ensure the proj4 string can fit in the HDF5 data type
    if( _srs_wkt.length() >= 256)
    {
        CHM_THROW_EXCEPTION(config_error,
                "Proj4 string needs to be < 256. Length: " + std::to_string(_srs_wkt.length()));
    }
    {  // Write the proj4
      H5::DataSpace dataspace(1, &proj4_dims);
      H5::Attribute attribute = file.createAttribute("/mesh/proj4", proj4_t, dataspace);
      attribute.write(proj4_t, _srs_wkt);
    }

    {  // Write the mesh version
        H5::DataSpace dataspace(1, &meshver_dims);
        H5::Attribute attribute = file.createAttribute("/mesh/version", meshver_t, dataspace);
        attribute.write(meshver_t, _version.to_string());
    }

    { // Write the partition method
        H5::DataSpace dataspace(1, &partition_type_dims);
        H5::Attribute attribute = file.createAttribute("/mesh/partition_method", partition_type_t, dataspace);
        attribute.write(partition_type_t, _partition_method);
    }

    {  // Write the is_geographic
      H5::DataSpace dataspace(1, &geographic_dims);
      H5::Attribute attribute = file.createAttribute("/mesh/is_geographic", PredType::NATIVE_HBOOL, dataspace);
      attribute.write(PredType::NATIVE_HBOOL, &_is_geographic);
    }

    {
        bool _is_partition = false; // will only ever be false here, true when created with the hd5 code path in partition
        H5::DataSpace dataspace(1, &partition_dims);
        H5::Attribute attribute = file.createAttribute("/mesh/is_partition", PredType::NATIVE_HBOOL, dataspace);
        attribute.write(PredType::NATIVE_HBOOL, &_is_partition);
    }

  } // end of try block



    // catch failure caused by the H5File operations
  catch (FileIException& error) {
    error.printErrorStack();
  }

  // catch failure caused by the DataSet operations
  catch (DataSetIException& error) {
    error.printErrorStack();
  }

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Parameters
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

  std::string par_filename = filename_base + "_param.h5";

  try {
    // Turn off the auto-printing when failure occurs so that we can
    // handle the errors appropriately
    Exception::dontPrint();

    H5::H5File file(par_filename,H5F_ACC_TRUNC);
    H5::Group group(file.createGroup("/parameters"));

    hsize_t ntri= size_global_faces();

    for (auto &par_iter : _parameters) {

      H5::DataSpace dataspace(1, &ntri);
      auto globalIDs = get_global_IDs();
      std::string par_location = "/parameters/" + par_iter;
      H5::DataSet dataset = file.createDataSet(par_location, PredType::NATIVE_DOUBLE, dataspace);

      std::vector<double> values(ntri);
#pragma omp parallel for
      for(size_t i=0; i<ntri; ++i) {
	auto face = _faces.at(i);
	values[i] = face->parameter(par_iter);
      }
      dataset.write(values.data(), PredType::NATIVE_DOUBLE);
    }

  } // end try block

    // catch failure caused by the H5File operations
  catch (FileIException& error) {
    error.printErrorStack();
  }

  // catch failure caused by the DataSet operations
  catch (DataSetIException& error) {
    error.printErrorStack();
  }

}

void attr_op(H5::H5Location &loc, const std::string attr_name,
             void *operator_data) {
  std::cout << attr_name << std::endl;
}

/*
 * Operator function.
 */
typedef struct _MeshParameters {
std::vector<std::string> names;
} MeshParameters;

herr_t
group_info(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata)
{
    hid_t group;
    MeshParameters *pars = (MeshParameters*)opdata;
    std::string str(name);
    group = H5Gopen2(loc_id, name, H5P_DEFAULT);
    // cout << "Name : " << str << endl; // Display the group name.
    (pars->names).push_back(str);
    H5Gclose(group);
    return 0;
}

void triangulation::load_mesh_from_h5(const std::string& mesh_filename)
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
        catch(AttributeIException& e)
        {
            v = ""; // will cause a default to version 1.0.0
        }
        _version.from_string(v);

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

    { // read the is_geographic
        H5::DataSpace dataspace(1, &geographic_dims);
        H5::Attribute attribute = file.openAttribute("/mesh/is_geographic");
        attribute.read(PredType::NATIVE_HBOOL, &_is_geographic);
    }

    { // are we loading from a h5 file that is participating in a partitioned mesh?
        try
        {

            H5::DataSpace dataspace(1, &partition_dims);
            H5::Attribute attribute = file.openAttribute("/mesh/is_partition");
            attribute.read(PredType::NATIVE_HBOOL, &_mesh_is_from_partition);


        }
        catch (AttributeIException& e)
        {
            // non partition meshes won't have this
            _mesh_is_from_partition = false;
        }

        SPDLOG_DEBUG("Loaded mesh is partitioned = {}",_mesh_is_from_partition);
    }

    if(!_mesh_is_from_partition && !_version.mesh_ver_meets_min_h5())
    {
        CHM_THROW_EXCEPTION(mesh_error, "h5 mesh version too old");
    }

    if(_mesh_is_from_partition && !_version.mesh_ver_meets_min_partition())
    {
        CHM_THROW_EXCEPTION(mesh_error, "partitioned mesh version too old");
    }

    try
    {

        {
            // Read the parition method
            H5::DataSpace dataspace(1, &partition_type_dims);
            H5::Attribute attribute = file.openAttribute("/mesh/partition_method");
            attribute.read(partition_type_t, _partition_method);

        }

        // if we are from a partition we are loading _local_sizes from the partition file not from the h5
        // they're almost always identical /unless/ we have created a subset debugging partition file (e.g., -v 1 -v 2)
        // in this case, the partition file will have the "real" number of partitions whereas this h5 is
        // before we subset it
        if(!_mesh_is_from_partition)
        {
            H5::DataSet dataset = file.openDataSet("/mesh/local_sizes");
            H5::DataSpace dataspace = dataset.getSpace();
            hsize_t nelem;
            int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);
            _local_sizes.resize(nelem);
            dataset.read(_local_sizes.data(), PredType::NATIVE_INT);
        }
    } catch(...)
    {
        CHM_THROW_EXCEPTION(mesh_error, "This h5 was produced by an older version of partition/meshpermutation.py and is lacking a key field. Please rerun these tools");
    }

    // CHM now needs pre-partitioning via metis method, any other partitioned methods are no longer supported
    if (_mesh_is_from_partition && _partition_method != "metis")
    {
        CHM_THROW_EXCEPTION(mesh_error, "CHM requires partitioned meshes built using the metis method");
    }

    {
        H5::DataSet dataset = file.openDataSet("/mesh/cell_global_id");
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t nelem;
        int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);
        _global_IDs.resize(nelem);
        dataset.read(_global_IDs.data(), PredType::NATIVE_INT);
    }

    std::vector<int> owner; //what MPIrank owns each triangle,
    {
        // if we have a h5 that isn't partitioned (ie mesh <v3.0.0) then we are entirely owned by rank 0
        if(_version.to_string() == "1.0.0" ||
            _version.to_string() == "1.2.0" ||
            _version.to_string() == "2.0.0")
        {
            owner.resize(_global_IDs.size(), 0);
        }
        else
        {
            H5::DataSet dataset = file.openDataSet("/mesh/owner");
            H5::DataSpace dataspace = dataset.getSpace();
            hsize_t nelem;
            int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);
            owner.resize(nelem);
            dataset.read(owner.data(), PredType::NATIVE_INT);
        }

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

            Point_3 pt(vertex[i][0], vertex[i][1], vertex[i][2]); // x y z
            _max_z = std::max(_max_z, vertex[i][2]);
            _min_z = std::min(_min_z, vertex[i][2]);

            _bounding_box.x_max = std::max(_bounding_box.x_max, vertex[i][0]);
            _bounding_box.x_min = std::min(_bounding_box.x_min, vertex[i][0]);

            _bounding_box.y_max = std::max(_bounding_box.y_max, vertex[i][1]);
            _bounding_box.y_min = std::min(_bounding_box.y_min, vertex[i][1]);

            Vertex_handle Vh = this->create_vertex();
            Vh->set_point(pt);
            Vh->set_id(i);
            _vertexes.push_back(Vh);

            // std::cout << i << ":";
            // for(int j=0; j< 3; ++j){
            //   std::cout << " " << vertex[i][j];
            // }
            // std::cout << "\n";
        }

        _num_vertex = _vertexes.size();
        SPDLOG_DEBUG("# nodes created = {}",_num_vertex);
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
            // get the global ID from file, so-as to support either pre partitioned or non partitioned meshes
            face->cell_global_id = _global_IDs.at(i);
            face->owner = owner.at(i);

            // this local id will include local ids for ghosts. However, that will need to be reset once partition
            // splits out the ghosts
            face->cell_local_id = i;

            // if we load from a partition, we need a way of linking the non contiguous global ids with the local ids
            if(_mesh_is_from_partition)
            {
                _global_to_locally_owned_index_map[_global_IDs.at(i)] = i;
            }

            if (_is_geographic)
            {
                face->_is_geographic = true;
            }

            face->_debug_ID =
                -(i +
                  1); // all ids will be negative starting at -1. Named ids (for output) will be positive starting at 0
            face->_debug_name = std::to_string(i);
            face->_domain = this;

            vert1->set_face(face);
            vert2->set_face(face);
            vert3->set_face(face);

            _faces.push_back(face);

        }
        _num_faces = _faces.size();
        SPDLOG_DEBUG("Created a mesh with {} triangles", this->size_faces());
    }

    {

        DataSet dataset = file.openDataSet("/mesh/neighbor");
        DataSpace dataspace = dataset.getSpace();

        hsize_t nelem;
        int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);

        if( elem.size() != nelem)
        {
            CHM_THROW_EXCEPTION(config_error,
                "Expected: " + std::to_string(elem.size()) + " neighborlists, got: " + std::to_string(nelem));
        }

        // Ensure enough space in the vector
        neigh.resize(nelem);
        // Default args read all of the dataspace
        dataset.read(neigh.data(), neighbor_t);

        SPDLOG_DEBUG("Building face neighbors");

        for (size_t i=0;i<nelem;i++)
        {

            auto face = _faces.at(i);

            // LOG_DEBUG << neigh[i][0] << " " << neigh[i][1] << " " << neigh[i][2];

            int neigh_i_0_idx = neigh[i][0];
            int neigh_i_1_idx = neigh[i][1];
            int neigh_i_2_idx = neigh[i][2];

//            // if we are loading from a partition, convert these to local ids
            if(_mesh_is_from_partition)
            {
                // neighbours that are missing come in as -1 so they won't cleanly remap

                neigh_i_0_idx = neigh_i_0_idx != -1 ? _global_to_locally_owned_index_map[neigh_i_0_idx] : -1;
                neigh_i_1_idx = neigh_i_1_idx != -1 ? _global_to_locally_owned_index_map[neigh_i_1_idx] : -1;
                neigh_i_2_idx = neigh_i_2_idx != -1 ? _global_to_locally_owned_index_map[neigh_i_2_idx] : -1;
            }

            if(    neigh_i_0_idx > static_cast<int>(nelem)
                   || neigh_i_1_idx > static_cast<int>(nelem)
                   || neigh_i_2_idx > static_cast<int>(nelem))
            {
                CHM_THROW_EXCEPTION(config_error, "Face " + std::to_string(i) + " has out of bound neighbors.");
            }

            //-1 is now the no neighbor value
            Face_handle face0 =  neigh.at(i)[0] != -1 ?_faces.at( neigh_i_0_idx) : nullptr;
            Face_handle face1 =  neigh.at(i)[1] != -1 ?_faces.at( neigh_i_1_idx) : nullptr;
            Face_handle face2 =  neigh.at(i)[2] != -1 ?_faces.at( neigh_i_2_idx) : nullptr;

            if( face0 != nullptr && face0->cell_global_id == face->cell_global_id)
            {
                SPDLOG_ERROR("At idx={} Face0 is trying to set itself as neighbour!", i);
                CHM_THROW_EXCEPTION(config_error, "Face0 is face!");
            }

            if( face1 != nullptr && face1->cell_global_id == face->cell_global_id)
            {
                SPDLOG_ERROR("At idx={} Face1 is trying to set itself as neighbour!", i);
                CHM_THROW_EXCEPTION(config_error, "Face1 is face!");
            }

            if( face2 != nullptr && face2->cell_global_id == face->cell_global_id)
            {
                SPDLOG_ERROR("At idx={} Face2 is trying to set itself as neighbour!", i);
                CHM_THROW_EXCEPTION(config_error, "Face2 is face!");
            }

            face->set_neighbors(face0,face1,face2);

        }
    }
}
void triangulation::_build_dDtree()
{
    SPDLOG_DEBUG("Building dD tree");

    size_t nfaces =   _faces.size();

    std::vector<Point_2> center_points(nfaces);

#pragma omp parallel for
    for(size_t ii=0; ii < nfaces; ++ii)
    {
        auto face = _faces.at(ii);
        Point_2 pt2(face->center().x(),face->center().y());
        center_points[ii] = pt2;

    }

    //make the search tree
    dD_tree = boost::make_shared<Tree>(boost::make_zip_iterator(boost::make_tuple( center_points.begin(),_faces.begin() )),
                                       boost::make_zip_iterator(boost::make_tuple( center_points.end(),  _faces.end() ) )
    );
}

void triangulation::from_partitioned_hdf5(const std::string& partition_filename,
                                          bool only_load_params,
                                          boost::filesystem::path cwd
                                          )
{
    // we are loading a partition mesh
    auto partition = read_json(partition_filename);

    auto nranks = partition.get<size_t>("ranks");

    //this will also be set when we load the h5 and is_parition attribute is read
    // however set it here just in case
    _mesh_is_from_partition = true;

#ifdef USE_MPI
    if (nranks != _comm_world.size())
    {
         CHM_THROW_EXCEPTION(config_error, "The partitioned mesh was configured for nranks="+std::to_string(nranks)+ " but CHM was run with " + std::to_string(_comm_world.size()) + " ranks." );
    }
#endif

    auto max_ghost_distance = partition.get<double>("max_ghost_distance");
    _num_global_faces = partition.get<size_t>("num_global_faces");

    // local sizes of partitions
    std::vector<int> mesh_partition_sizes;
    try
    {
        for(auto &itr : partition.get_child("local_sizes"))
        {
	  auto size = stoi(itr.second.data());
	  mesh_partition_sizes.push_back(size);
        }
    }
    catch(pt::ptree_bad_path &e)
    {
        CHM_THROW_EXCEPTION(config_error, "A mesh is required and was not found");
    }
    _local_sizes.resize(nranks);
    for( size_t i=0; i < nranks; ++i)
    {
      _local_sizes.at(i) = mesh_partition_sizes.at(i);
    }



    // Paths for parameter files
    std::vector<std::string> mesh_file_paths;
    try
    {
        for(auto &itr : partition.get_child("meshes"))
        {
            auto path = boost::filesystem::path(itr.second.data());
            if (path.is_relative())
                mesh_file_paths.push_back( (cwd / path).string() );
            else
                mesh_file_paths.push_back(path.string() );
        }
    }
    catch(pt::ptree_bad_path &e)
    {
        CHM_THROW_EXCEPTION(config_error, "A mesh is required and was not found");
    }

    // Paths for parameter files
    std::vector<std::string> param_file_paths;
    try
    {
        for(auto &itr : partition.get_child("parameters"))
        {
            auto param_file = boost::filesystem::path(itr.second.data());
            if(param_file.is_relative())
                param_file_paths.push_back( (cwd /param_file).string() );
            else
                param_file_paths.push_back( param_file.string() );
        }
    }
    catch(pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG("No addtional parameters found in mesh section.");
    }

    // Paths for initial condition files
    std::vector<std::string> initial_condition_file_paths;
    try
    {
        for(auto &itr : partition.get_child("initial_conditions"))
        {
            auto ic_file = itr.second.data();
            SPDLOG_DEBUG("Found initial condition file: {}",ic_file);
            CHM_THROW_EXCEPTION(config_error,"Initial conditions not currently supported with partition tool.");
        }
    }
    catch(pt::ptree_bad_path &e)
    {
        SPDLOG_DEBUG("No addtional initial conditions found in mesh section.");
    }

    // keep just our rank's. The hdf5 loader expects a vector of potential files so make a vector of 1

    int tmp_rank = 0; // play nicely in non mpi builds
#ifdef USE_MPI
    tmp_rank = _comm_world.rank();
#endif

    auto mesh_file = mesh_file_paths[tmp_rank];

    std::vector<std::string> param_file;

    if(param_file_paths.size() > 0)
        param_file.push_back( param_file_paths[tmp_rank] );

    if(!only_load_params)
        from_hdf5(mesh_file, {}, {}, true );

    if(only_load_params)
        load_hdf5_parameters(param_file);

}
void triangulation::from_hdf5(const std::string& mesh_filename,
			      const std::vector<std::string>& param_filenames,
			      const std::vector<std::string>& ic_filenames,
                              bool delay_param_ic_load
			      )
{
    try
    {
        load_mesh_from_h5(mesh_filename);
    }
    // catch failure caused by the H5File operations
    catch (FileIException& e)
    {
        e.printErrorStack();
        CHM_THROW_EXCEPTION(mesh_error, "Error loading HDF5 file: " + mesh_filename);
    }
    // catch failure caused by the DataSet operations
    catch (DataSetIException& e) {
        e.printErrorStack();
    }

    // if we are partitioned, just load it from file
    if(_mesh_is_from_partition)
    {
        load_partition_from_mesh(mesh_filename);
    }
    else
    {
        // otherwise, compute it
        partition_mesh();

#ifdef USE_MPI
        determine_local_boundary_faces();
        determine_process_ghost_faces_nearest_neighbors();

        // TODO: Need to auto-determine how far to look based on module setups
        determine_process_ghost_faces_by_distance(100.0);
#endif
    }

#ifdef USE_MPI
    determine_ghost_owners();
    setup_nearest_neighbor_communication();

//    print_ghost_neighbor_info();

#endif // USE_MPI

    _build_dDtree();

    // load param
    if(!delay_param_ic_load)
        load_hdf5_parameters(param_filenames);

    // TODO: include initial condition files

}

void triangulation::load_hdf5_parameters( const std::vector<std::string>& param_filenames)
{
    for (auto param_filename : param_filenames)
    {
        try
        {
            // Turn off the auto-printing when failure occurs so that we can
            // handle the errors appropriately
            Exception::dontPrint();

            // Open an existing file and dataset.
            H5File file(param_filename, H5F_ACC_RDONLY);

            // Space for extracting the parameter info
            std::unique_ptr<MeshParameters> pars(new MeshParameters);

            // Extract all of the parameter info (names) from the file
            Group group = file.openGroup("parameters");
            herr_t idx = H5Literate(group.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, group_info, (void*)pars.get());
            for (auto const& name : pars->names)
            {
                _parameters.insert(name);
                // std::cout << "Here: " << name << "\n";
            }

            if(_mesh_is_from_partition)
            {
// init the parameter storage on each face
// as we are using a pre-partitioned mesh, _faces holds local+ghosts, so can do it in one go which is
// faster
#pragma omp parallel for
                for (size_t i = 0; i < _faces.size(); i++)
                {
                    _faces.at(i)->init_parameters(_parameters);
                }
            }
            else
            {
// if we are not reading from a partitioned file, we need to ensure we do the local faces + ghosts
// separetely
// init the parameter storage on each face
#pragma omp parallel for
                for (size_t i = 0; i < _num_faces; i++)
                {
                    face(i)->init_parameters(_parameters);
                }

// init the parameter storage for the ghost regions
#pragma omp parallel for
                for (size_t i = 0; i < _ghost_faces.size(); i++)
                {
                    _ghost_faces.at(i)->init_parameters(_parameters);
                }
            }
            // Data buffer for reading from file (before packing into faces)
            std::vector<double> data(_mesh_is_from_partition ? _faces.size() : _num_faces);

            // Read all of the parameters from file and store them in the faces
            for (auto const& name : pars->names)
            {

                SPDLOG_DEBUG("Applying parameter: {}",name);

                // _num_faces will have been set to the local face size in MPI mode
                // however, if we are reading from a partitioned mesh file, we can actually load the entire set of
                // faces' data in one go
                hsize_t local_size = static_cast<hsize_t>(
                    _mesh_is_from_partition ? _faces.size() : _num_faces);

                // if we are partitioned, then every file needs to start at offset 0
                hsize_t offset = static_cast<hsize_t>(_mesh_is_from_partition ? 0 : global_cell_start_idx);


                DataSet dataset = group.openDataSet(name);
                DataSpace dataspace = dataset.getSpace();

                hsize_t nelem;
                int ndims = dataspace.getSimpleExtentDims(&nelem);

                dataspace.selectHyperslab(H5S_SELECT_SET, &local_size, &offset);

                DataSpace memspace(1, &nelem);
                hsize_t zero = 0;
                memspace.selectHyperslab(H5S_SELECT_SET, &local_size, &zero);

                dataset.read(data.data(), PredType::NATIVE_DOUBLE, memspace, dataspace);

                if(_mesh_is_from_partition)
                {
// since the params are for all our faces + ghosts, we can load it all at once
#pragma omp parallel for
                    for (size_t i = 0; i < _faces.size(); i++)
                    {
                        _faces.at(i)->parameter(name) = data.at(i);
                    }
                }
                else
                {
// we have to load the ghosts and local faces separate.
#pragma omp parallel for
                    for (size_t i = 0; i < _num_faces; i++)
                    {
                        face(i)->parameter(name) = data[i];
                    }

                    SPDLOG_DEBUG(" Applying {} for ghost regions {} elements):",
                                  name,  _ghost_faces.size(), name);

                    // Read parameters for each ghost face individually
                    hsize_t one = 1;
                    // Do NOT do this loop in parallel (internal state of HDF5)
                    for (size_t i = 0; i < _ghost_faces.size(); i++)
                    {
                        auto face = _ghost_faces.at(i);
                        // when we are reading from a paritioned file, the ghosts will be below the non ghosts
                        hsize_t global_id = _mesh_is_from_partition ? _local_faces.size() + i : face->cell_global_id;

                        // Position and size in file
                        dataspace.selectHyperslab(H5S_SELECT_SET, &one, &global_id);
                        // Position and size in variable
                        memspace.selectHyperslab(H5S_SELECT_SET, &one, &zero);

                        double value;
                        dataset.read(&value, PredType::NATIVE_DOUBLE, memspace, dataspace);
                        face->parameter(name) = value;
                    }
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

    } // end of param_filenames loop

}

void triangulation::reorder_faces(std::vector<size_t> permutation)
{
  // NOTE: be careful with evaluating this, the 'cell_global_id's and a
  // cell's position in the '_faces' vec are unrelated. They must both
  // be modified (ie. renumber the 'cell_global_id's, AND sort the '_faces'
  // vector) before the new ordering is consistent.

  assert( permutation.size() == size_faces() );
  SPDLOG_DEBUG("Reordering faces");

  // Update the IDs on all faces
  #pragma omp parallel for
  for (size_t ind = 0; ind < permutation.size(); ++ind)
  {

	     size_t old_ID = permutation.at(ind);
	     size_t new_ID = ind;

	     auto face = _faces.at(old_ID);
	     face->cell_global_id = new_ID;

  }

  // Sort the faces in the new ordering
  tbb::parallel_sort(_faces.begin(), _faces.end(),
  		     [](triangulation::Face_handle fa, triangulation::Face_handle fb)->bool
  		     {
  		       return fa->cell_global_id < fb->cell_global_id;
  		     });
}

void triangulation::load_partition_from_mesh(const std::string& mesh_filename)
{
    // This differs from partition_mesh() in that the partitioned mesh has ghosts mixed in with the faces so that
    // the full geometry could be written to file. Like the rest of the partitioned this does not support a non MPI code path

    // Ensure that the local boundary faces have been determined, but ghost neighbors have not been set
    assert( _ghost_faces.size() == 0 );

#ifdef USE_MPI
    int my_rank = _comm_world.rank();

    SPDLOG_DEBUG("Loading partition");

    // Turn off the auto-printing when failure occurs so that we can
    // handle the errors appropriately
    Exception::dontPrint();
    // Open an existing file and dataset.
    H5File file(mesh_filename, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet("/mesh/ghost_type");
    H5::DataSpace dataspace = dataset.getSpace();
    hsize_t nelem;
    int ndims = dataspace.getSimpleExtentDims(&nelem, NULL);

    std::vector<int> ghost_info(nelem);
    dataset.read(ghost_info.data(), PredType::NATIVE_INT);
    file.close();

    SPDLOG_DEBUG(mesh_filename);

    // when we get to here, we have a mix of ghosts and not ghosts in _local_faces and don't have the sizes
    // we need to pick this apart
    // since the faces are sorted coming from the partition tool, they will be sorted here. Since we grow arrays we
    // can't do this in parallel which will preserve the sort. If this is every changed, then need to sort


    // here we loop through all (incl ghosts!) to figure out where everything should go.
    // DO NOT do this in parallel (at the moment) as it's not thread safe
    size_t local_face_i = 0;
    for (size_t i = 0; i < _faces.size(); ++i)
    {
        _faces[i]->ghost_type = ghost_info[i];

        // owner has already been set from loading in the h5
//        _faces[i]->owner = my_rank;
        if (ghost_info[i] != GHOST_TYPE::NONE)
        {
            _faces[i]->is_ghost = true;

            _ghost_faces.push_back(_faces[i]);
            if(ghost_info[i] == GHOST_TYPE::NEIGH)
                _ghost_neighbors.push_back(_faces[i]);

        } else
        {
            _faces[i]->is_ghost = false;

            // Although we set thsin from_h5 we need to reset it here such that it doesn't include ghosts
            _faces[i]->cell_local_id = local_face_i;

            _local_faces.push_back(_faces[i]);

            _global_to_local_faces_index_map[_faces[i]->cell_global_id] = local_face_i;
            local_face_i++;

        }
    }

    // Set up so that all processors know how 'big' all other processors are
    _num_faces_in_partition.resize(_comm_world.size());
    for (unsigned int i = 0; i < _comm_world.size(); ++i)
    {
        _num_faces_in_partition.at(i) = _local_sizes.at(i);
    }

    if( _local_faces.size() != _num_faces_in_partition.at(_comm_world.rank()) )
    {
        SPDLOG_ERROR("Expected to see #faces= {}\nbut got #faces={}",
                      _num_faces_in_partition.at(_comm_world.rank()), _local_faces.size());
        CHM_THROW_EXCEPTION(mesh_error, "local_size and partition size mismatch" );
    }

//    // each processor only knows its own start and end indices
//    size_t face_start_idx = 0;
//    size_t face_end_idx = _num_faces_in_partition[0] - 1;
//    for (int i = 1; i <= _comm_world.rank(); ++i)
//    {
//        face_start_idx += _num_faces_in_partition[i - 1];
//        face_end_idx += _num_faces_in_partition[i];
//    }
//
//    if (face_start_idx != _local_faces.front()->cell_global_id ||
//        face_end_idx != _local_faces.back()->cell_global_id)
//    {
//        LOG_ERROR << "The computed and read start/end index don't match. Computed:\n"
//                  << "\tface_start_idx = " << face_start_idx << "\n"
//                  << "\tface_end_idx = " << face_end_idx << "\n"
//                  << "Read:\n"
//                  << "\t_local_faces[0] = " << _local_faces.front()->cell_global_id << "\n"
//                  << "\t_local_faces[-1] = " << _local_faces.back()->cell_global_id;
//        CHM_THROW_EXCEPTION(mesh_error, "Computed and read global indexs don't match");
//    }

    // Store in the public members
    global_cell_start_idx = _local_faces.front()->cell_global_id;
    global_cell_end_idx   = _local_faces.back()->cell_global_id;

    _num_faces = _local_faces.size();

    // _global_IDs must contain (in the same order) cell_global_id for the faces
    // in _local_faces
    {
      std::vector<int> tmp(_local_faces.size());
      std::transform(_local_faces.begin(), _local_faces.end(),
		     tmp.begin(),
		     [](mesh_elem e) { return e->cell_global_id; });
      std::swap(tmp, _global_IDs);
    }


    SPDLOG_DEBUG( "MPI Process {}: start {}, end {} , number {}",
                  my_rank, global_cell_start_idx, global_cell_end_idx, _local_faces.size());

#endif
}

void triangulation::partition_mesh_nonMPI(size_t _num_global_faces)
{
    // If we are not using MPI, some code paths might still want to make use of these
    //  initialized to [0, num_faces - 1]
    global_cell_start_idx = 0;
    global_cell_end_idx = _num_faces - 1;

    _global_IDs.resize(_num_global_faces);
#pragma omp parallel for
    for (size_t i = 0; i < _num_global_faces; ++i)
    {
        _global_IDs[i] = i;
        _faces.at(i)->is_ghost = false;
        _faces.at(i)->cell_local_id =
            i; // Mesh has been (potentially) reordered before this point. Set the local_id correctly
    }

    // make sure these setup when in MPI mode for partition
    _num_faces =  _faces.size();
    _num_global_faces = _faces.size();
    _local_faces = _faces;

    SPDLOG_DEBUG("Face numbering : start 0, end {}, number {}",(_num_global_faces - 1), _local_faces.size());
}

void triangulation::partition_mesh()
{
    /*
      Basically, this first attempt splits the number of mesh points as equally as possible
    */
    // 1. Determine number of (processor) locally owned faces
    // 2. Determine (processor) locally owned indices of faces
    // 3. Set locally owned is_ghost=false

    _num_global_faces = _faces.size();

#ifdef USE_MPI

    SPDLOG_DEBUG("Partitioning mesh");

    int my_rank = _comm_world.rank();

    // Set up so that all processors know how 'big' all other processors are
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

    // Loop can't be parallel due to modifying map
    for (size_t local_ind = 0; local_ind < _local_faces.size(); ++local_ind)
    {
        size_t offset_idx = global_cell_start_idx + local_ind;
        _global_to_locally_owned_index_map[_global_IDs.at(offset_idx)] = local_ind;

        _faces.at(_global_IDs.at(offset_idx))->is_ghost = false;
        _faces.at(_global_IDs.at(offset_idx))->owner = my_rank;

        //this was set in from_h5 but needs to be reset to not include any of the ghosts
        _faces.at(_global_IDs.at(offset_idx))->cell_local_id = local_ind;
        _local_faces.at(local_ind) = _faces.at(_global_IDs.at(offset_idx));
    }

    // _global_IDs must contain (in the same order) cell_global_id for the faces
    // in _local_faces
    {
      std::vector<int> tmp(_local_faces.size());
      std::transform(_local_faces.begin(), _local_faces.end(),
		     tmp.begin(),
		     [](mesh_elem e) { return e->cell_global_id; });
      std::swap(tmp, _global_IDs);
    }

    _num_faces = _local_faces.size();
    SPDLOG_DEBUG( "MPI Process {}: start {}, end  {} , number {}",
                  my_rank, global_cell_start_idx,  global_cell_end_idx, _local_faces.size());

#else // do not USE_MPI

    partition_mesh_nonMPI(_num_global_faces);

#endif // USE_MPI
}

void triangulation::determine_local_boundary_faces()
{
  /*
    - Store handles to boundary faces on the locally owned process
    - Also store boolean value "is_global_boundary"
  */

  using th_safe_multicontainer_type = std::vector< std::pair<mesh_elem,bool> >[];

#ifdef USE_MPI
  // Need to ensure we're starting from nothing?
  assert( _boundary_faces.size() == 0 );

  SPDLOG_DEBUG("Determining local boundary faces");

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
  SPDLOG_DEBUG("MPI Process {}  has {} boundary faces.", _comm_world.rank(), _boundary_faces.size());

  SPDLOG_DEBUG("MPI Process {}  _faces.size():{} _local_faces.size():{}", _comm_world.rank(), _faces.size(), _local_faces.size());

  // _comm_world.barrier();
  // exit(0);

#endif
}

void triangulation::determine_process_ghost_faces_nearest_neighbors()
{
  // NOTE that this algorithm is not implemented for multithread
  // - multithread can be implemented similarly to determine_local_boundary_faces
  //   - not a priority since number of boundary faces should be small

  // Ensure that the local boundary faces have been determined, but ghost nearest neighbors have not been set
  assert( _boundary_faces.size() != 0 );
  assert( _ghost_neighbors.size() == 0 );

  SPDLOG_DEBUG("Determining ghost region info");

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

#ifdef USE_MPI
  SPDLOG_DEBUG("MPI Process {} has {} ghosted nearest neighbors.", _comm_world.rank(), _ghost_neighbors.size());
#endif



}

// Generate unique tags for isend/irecv
// - unique id for proc m communicating with proc p
// - tags will be unique for up to 9999 MPI processes
// - message between m and p must be completed before a new one started
// - this is ensured with mpi::wait_all, as in ghost_neighbors_communicate_variable
int generate_unique_send_tag(int my_rank, int partner_rank){
  return 10000*my_rank + partner_rank;
}
int generate_unique_recv_tag(int my_rank, int partner_rank){
  return my_rank + 10000*partner_rank;
}

void triangulation::determine_ghost_owners()
{
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
        _ghost_neighbor_owners[i] = _ghost_neighbors.at(i)->owner;

//        SPDLOG_DEBUG("Inc global ind {}",global_ind);
//        _ghost_neighbor_owners[i] = determine_owner_of_global_index(global_ind,
//                                                                    _num_faces_in_partition);
        // local faces are set to current mpirank, so set the ghosts to be ID"d who owns them
//        _ghost_neighbors[i]->owner = _ghost_neighbor_owners[i];

        // on first it, no value of prev_owner exists... set it
        if(i==0) prev_owner = _ghost_neighbor_owners[i];

        // if owner different from last owner, store prev segment's ownership info
        if (prev_owner != _ghost_neighbor_owners[i])
        {
            num_partners++;
            _comm_partner_ownership[prev_owner] = std::make_pair(start_index, i-start_index);
            start_index=i;
        } else if (i ==_ghost_neighbors.size()-1) {
            _comm_partner_ownership[prev_owner] = std::make_pair(start_index, i-start_index+1);
        }

	// If the last neighbor rank owns only a single ghost
	// - would have followed the if branch of previous conditional
	//   - this implies i==start_index
	// - we set the last ownership info explicitly
	if (start_index ==_ghost_neighbors.size()-1) {
	  assert(i==start_index);
	  _comm_partner_ownership[_ghost_neighbor_owners[i]] = std::make_pair(i,1);
	}

        // prep prev_owner for next iteration
        prev_owner=_ghost_neighbor_owners[i];
    }

    for(size_t i=0; i<_ghost_neighbors.size(); ++i)
    {
        _global_index_to_local_ghost_map[_ghost_neighbors[i]->cell_global_id] = static_cast<int>(i);
    }

#ifdef USE_MPI
    SPDLOG_DEBUG("MPI Process {}  has {} communication partners:", _comm_world.rank(), _comm_partner_ownership.size());
    for (auto it : _comm_partner_ownership)
    {
        SPDLOG_DEBUG("  # Process {} partner {} owns local ghost indices {} to {}",  _comm_world.rank(), it.first, it.second.first, it.second.first+it.second.second-1);
    }

#endif
}
void triangulation::setup_nearest_neighbor_communication()
{

// Function is meaningful only when using MPI
#ifdef USE_MPI

  /*
    Each process knows what global IDs it needs, and which other process owns them
    - need to let those other processes know which indices to send to us
   */

  std::vector<boost::mpi::request> reqs;

  // Send the list of needed global indices to the required processes
  for( auto it : _comm_partner_ownership ) {

    auto partner_id = it.first;
    // LOG_DEBUG << "MPI Process " << _comm_world.rank() << " partner_id: " << partner_id;
    auto start_idx = it.second.first;
    auto length = it.second.second;

    // Iterate over global indices
    std::vector<int> id_indices(_ghost_neighbors.size());
    std::transform(_ghost_neighbors.begin(),_ghost_neighbors.end(),id_indices.begin(),
		   [](mesh_elem e){ return e->cell_global_id; });
    // Copy subvector of indices to communicate
    ghost_indices_to_recv[partner_id].resize(length);
    std::copy(id_indices.begin()+start_idx,
    	      id_indices.begin()+start_idx+length,
              ghost_indices_to_recv[partner_id].begin());

    // for (auto it : id_indices)
    //   {
    // 	LOG_DEBUG << "  # Process " << _comm_world.rank() << " id_index " << it;
    //   }



    // Copy relevant portion of ghost neighbors to the "local_faces_to_recv"
    ghost_faces_to_recv[partner_id].resize(length);
    std::copy(_ghost_neighbors.begin()+start_idx,
	      _ghost_neighbors.begin()+start_idx+length,
              ghost_faces_to_recv[partner_id].begin());
    // for (auto it : sub_indices)
    //   {
    // 	LOG_DEBUG << "  # Process " << _comm_world.rank() << " partner " << partner_id << " global_id " << it;
    //   }

    // Send indices
    int send_tag = generate_unique_send_tag(_comm_world.rank(), partner_id);
    reqs.push_back(_comm_world.isend(partner_id, send_tag, ghost_indices_to_recv[partner_id]));

  }

  // Recv the list of global indices to send to the required processes
  for( auto it : _comm_partner_ownership ) {

    auto partner_id = it.first;
    // LOG_DEBUG << "MPI Process " << _comm_world.rank() << " partner_id: " << partner_id;

    // Recv indices
    // Note: opposite constants from send tags
    int recv_tag = generate_unique_recv_tag(_comm_world.rank(), partner_id);

    // Receive directly into the global_indices_to_send map
    reqs.push_back(_comm_world.irecv(partner_id, recv_tag, global_indices_to_send[partner_id] ));

  }

  // Wait for all comms to finish before proceeding
  boost::mpi::wait_all(reqs.begin(), reqs.end());

  for (auto it : global_indices_to_send) {
    auto partner_id = it.first;
    auto indices = it.second;
    // Get local indices that need to be sent
    local_indices_to_send[partner_id].resize(indices.size());
    std::transform(indices.begin(), indices.end(),
		   local_indices_to_send[partner_id].begin(),
		   [this](int ind){ return _global_to_local_faces_index_map[ind]; });
    // Get local faces that need to be sent
    local_faces_to_send[partner_id].resize(indices.size());
    std::transform(local_indices_to_send[partner_id].begin(), local_indices_to_send[partner_id].end(),
		   local_faces_to_send[partner_id].begin(),
		   [this](int ind) {return _local_faces.at(ind); });
  }

  // for (auto it : global_indices_to_send)
  //   {
  //     auto partner = it.first;
  //     for (auto it_glob : it.second) {
  //   	LOG_DEBUG << "  # Process " << _comm_world.rank() << " partner " << partner << " global_id to send " << it_glob;
  //     }
  //   }

  // _comm_world.barrier();
  // exit(0);

#endif // USE_MPI

}

void triangulation::print_ghost_neighbor_info()
{
  // Simple text file output, separate files for each MPI rank
#ifdef USE_MPI
  auto myrank = _comm_world.rank();
  std::string filename = "ghost_neighbor_info_rank_" + std::to_string(myrank);

  std::ofstream outfile(filename);

  SPDLOG_DEBUG("Rank {} writing ghost neighbor info to file.", myrank);

  outfile << "#(position in ghost array) (cell_global_id) (owner)\n";
  for (int ii=0; ii < _ghost_neighbors.size(); ++ii) {
    outfile << ii << " "
	    << _ghost_neighbors[ii]->cell_global_id << " "
	    << _ghost_neighbor_owners[ii] << "\n";
  }
  outfile.close();
#endif // USE_MPI
}

void triangulation::ghost_neighbors_communicate_variable(const std::string& var)
{
    // This supports the use case if _s no-oped to const char * via
    uint64_t hash = xxh64::hash (var.c_str(), var.length());
//    LOG_DEBUG << hash << " " << var;
    ghost_neighbors_communicate_variable(hash);
}

void triangulation::ghost_neighbors_communicate_variable(const uint64_t& var)
{

// Function is meaningful only when using MPI
#ifdef USE_MPI

  // For each communication partner:
  // - pack vectors of the variable to send
  // - send/recv it
  // - unpack the recv'd vectors into mesh_elem->face_data (so it can be used exactly as local info)

  std::vector<boost::mpi::request> reqs;

  for(auto it : local_faces_to_send) {
    auto partner_id = it.first;
    auto faces = it.second;
    std::vector<double> send_buffer(faces.size());
    std::transform(faces.begin(), faces.end(),
		   send_buffer.begin(),
		   [var](mesh_elem e){
		     return (*e)[var]; });

    for(int i=0; i < send_buffer.size(); ++i) {
      double val = send_buffer[i];
      if( isnan(val) )	{
	  auto f = local_faces_to_send[partner_id][i];
	  SPDLOG_DEBUG("-------------------------------------------------");
	  SPDLOG_DEBUG("Detected SEND variable is NaN:");
	  SPDLOG_DEBUG("\tmy rank:            {}",f->owner);
	  SPDLOG_DEBUG("\tdestination rank:  {}",partner_id);
	  SPDLOG_DEBUG("\tsend_buffer entry: {}",i);
	  SPDLOG_DEBUG("\tcell_global_id:    {}",f->cell_global_id);
	  SPDLOG_DEBUG("\tcell_local_id:      {}",f->cell_local_id);
          SPDLOG_DEBUG("\tvalue:       {}",val);
//	  _mpi_env.abort(-1);
	}
    }


    // Send variables
    int send_tag = generate_unique_send_tag(_comm_world.rank(), partner_id);
    _comm_world.isend(partner_id, send_tag, send_buffer);

  }

  // map of received data from comm partners
  std::map< int, std::vector<double>> recv_buffer;
  for( auto it : ghost_faces_to_recv) {
    auto partner_id = it.first;
    auto faces = it.second;
    recv_buffer[partner_id] = std::vector<double>(faces.size(),0.0);
  }

  for(auto it : ghost_faces_to_recv) {
    auto partner_id = it.first;

    // Note only recv gets added to the list of requests to watch for
    int recv_tag = generate_unique_recv_tag(_comm_world.rank(), partner_id);
    reqs.push_back(_comm_world.irecv(partner_id, recv_tag, recv_buffer[partner_id] ));
  }

  // Wait for all communication to me before proceeding
  boost::mpi::wait_all(reqs.begin(), reqs.end());

  // Check for NaN entries in recv_buffer
  for(auto it : recv_buffer) {
    auto partner_id = it.first;
    auto values = it.second;
    for(int i=0; i < values.size(); ++i) {
      double val = values[i];
      if( isnan(val) )	{
	auto f = ghost_faces_to_recv[partner_id][i];
	  SPDLOG_DEBUG("-------------------------------------------------");
	  SPDLOG_DEBUG("Detected RECV variable is NaN:");
	  SPDLOG_DEBUG("\tmy rank:            {}",f->owner);
	  SPDLOG_DEBUG("\tsent from rank:     {}",partner_id);
	  SPDLOG_DEBUG("\trecv_buffer entry:  {}",i);
	  SPDLOG_DEBUG("\tcell_global_id:     {}",f->cell_global_id);
	  SPDLOG_DEBUG("\tcell_local_id:       {}",f->cell_local_id);
          SPDLOG_DEBUG("\tvalue:       {}",val);
//	  _mpi_env.abort(-1);
	}
    }

  }

  // Pack the data into the face pointers
  for(auto it : ghost_faces_to_recv) {
    auto partner_id = it.first;
    auto faces = it.second;

    auto recv_it = recv_buffer[partner_id].begin();
    for (auto f : faces ) {

      (*f)[var] = *recv_it;
//      LOG_DEBUG << "Rank=" <<f->owner << " Var="<<var<<" *f="<<(*f)[var]<< " recv_it="<<*recv_it << " is_ghost="<<f->is_ghost;
      ++recv_it;
    }

  }

#endif // USE_MPI

}

void triangulation::ghost_to_neighbors_communicate_variable(const std::string& var)
{
    // This supports the use case if _s no-oped to const char * via
    uint64_t hash = xxh64::hash (var.c_str(), var.length());
    ghost_to_neighbors_communicate_variable(hash);
}

void triangulation::ghost_to_neighbors_communicate_variable(const uint64_t& var)
{

// Function is meaningful only when using MPI
#ifdef USE_MPI

    // For each communication partner:
    // - pack vectors of the variable to send
    // - send/recv it
    // - unpack the recv'd vectors into mesh_elem->face_data (so it can be used exactly as local info)

    std::vector<boost::mpi::request> reqs;

    for(auto it : ghost_faces_to_recv) {
        auto partner_id = it.first;
        auto faces = it.second;
        std::vector<double> send_buffer(faces.size());
        std::transform(faces.begin(), faces.end(),
                       send_buffer.begin(),
                       [var](mesh_elem e){
                           return (*e)[var]; });

        for(int i=0; i < send_buffer.size(); ++i) {
            double val = send_buffer[i];
            if( isnan(val) )	{
                auto f = ghost_faces_to_recv[partner_id][i];
                SPDLOG_DEBUG("-------------------------------------------------");
                SPDLOG_DEBUG("Detected SEND variable is NaN:");
                SPDLOG_DEBUG("\tmy rank:            {}",f->owner);
                SPDLOG_DEBUG("\tdestination rank:  {}",partner_id);
                SPDLOG_DEBUG("\tsend_buffer entry: {}",i);
                SPDLOG_DEBUG("\tcell_global_id:    {}",f->cell_global_id);
                SPDLOG_DEBUG("\tcell_local_id:      {}",f->cell_local_id);
                //	  _mpi_env.abort(-1);
            }
        }


        // Send variables
        int send_tag = generate_unique_send_tag(_comm_world.rank(), partner_id);
        _comm_world.isend(partner_id, send_tag, send_buffer);

    }

    // map of received data from comm partners
    std::map< int, std::vector<double>> recv_buffer;
    for( auto it : local_faces_to_send) {
        auto partner_id = it.first;
        auto faces = it.second;
        recv_buffer[partner_id] = std::vector<double>(faces.size(),0.0);
    }

    for(auto it : local_faces_to_send) {
        auto partner_id = it.first;

        // Note only recv gets added to the list of requests to watch for
        int recv_tag = generate_unique_recv_tag(_comm_world.rank(), partner_id);
        reqs.push_back(_comm_world.irecv(partner_id, recv_tag, recv_buffer[partner_id] ));
    }

    // Wait for all communication to me before proceeding
    boost::mpi::wait_all(reqs.begin(), reqs.end());

    // Check for NaN entries in recv_buffer
    for(auto it : recv_buffer) {
        auto partner_id = it.first;
        auto values = it.second;
        for(int i=0; i < values.size(); ++i) {
            double val = values[i];
            if( isnan(val) )	{
                auto f = local_faces_to_send[partner_id][i];
                SPDLOG_DEBUG("-------------------------------------------------");
                SPDLOG_DEBUG("Detected RECV variable is NaN:");
                SPDLOG_DEBUG("\tmy rank:            {}",f->owner);
                SPDLOG_DEBUG("\tsent from rank:     {}",partner_id);
                SPDLOG_DEBUG("\trecv_buffer entry:  {}",i);
                SPDLOG_DEBUG("\tcell_global_id:     {}",f->cell_global_id);
                SPDLOG_DEBUG("\tcell_local_id:       {}",f->cell_local_id);
                //	  _mpi_env.abort(-1);
            }
        }

    }

    // Pack the data into the face pointers
    for(auto it : local_faces_to_send) {
        auto partner_id = it.first;
        auto faces = it.second;

        auto recv_it = recv_buffer[partner_id].begin();
        for (auto f : faces ) {
            (*f)[var] = *recv_it;
            ++recv_it;
        }

    }

#endif // USE_MPI

}

void dfs_to_max_distance_aux(mesh_elem starting_face, double max_distance, mesh_elem face, std::unordered_set<mesh_elem> &visited)
{
  // DFS auxiliary function, does all of the work constructing the set of faces

  // if we're outside of the distance or already visited, we're done
  bool is_not_within_distance = (math::gis::distance(starting_face->center(),face->center()) > max_distance);
  bool is_visited             = (visited.find(face) != visited.end());
  if(  is_not_within_distance || is_visited  )
  {
    return;
  }

  // otherwise, visit face, and move on to neighbors
  visited.insert(face);
  for(int i = 0; i < 3; ++i)
  {
      auto neigh = face->neighbor(i);
      if(neigh != nullptr)
      {
	dfs_to_max_distance_aux(starting_face, max_distance, neigh, visited);
      }
  }
}

std::vector<mesh_elem> dfs_to_max_distance(mesh_elem starting_face, double max_distance)
{
  // Depth first search out to a maximum distance
  // - unknown number of graph elements

  // if a mesh_elem appears in this set, it has been visited
  std::unordered_set<mesh_elem> visited;
  dfs_to_max_distance_aux(starting_face, max_distance, starting_face, visited);
  // Convert to vector
  std::vector<mesh_elem> visited_vec(std::begin(visited), std::end(visited));
  return visited_vec;
}

void triangulation::determine_process_ghost_faces_by_distance(double max_distance)
{
  // NOTE that this algorithm is not implemented for multithread
  // - multithread can be implemented similarly to determine_local_boundary_faces
  //   - not a priority since number of boundary faces should be small

  // Ensure that the local boundary faces have been determined, but ghost neighbors have not been set
  assert( _boundary_faces.size() != 0 );
  assert( _ghost_faces.size() == 0 );

  // Vector for append speed
  std::vector< mesh_elem > ghosted_boundary_neighbors;

  for(size_t face_index=0; face_index< _boundary_faces.size(); ++face_index)
  {
    // face_index is a local index... get the face handle
    auto face = _boundary_faces.at(face_index).first;
    // separate the ghosted and non-ghosted neighbors
    // std::vector<mesh_elem> current_neighbors = find_faces_in_radius(face->center().x(),face->center().y(), max_distance);
    std::vector<mesh_elem> current_neighbors = dfs_to_max_distance(face, max_distance);
    auto pivot = std::partition(std::begin(current_neighbors),std::end(current_neighbors),
    				[] (mesh_elem neigh) {
    				  return neigh->is_ghost == true;
    				});
    current_neighbors.erase(pivot,std::end(current_neighbors));

    ghosted_boundary_neighbors.insert(std::end(ghosted_boundary_neighbors),
    				       current_neighbors.begin(),current_neighbors.end());
  }

  // Convert to a set to remove duplicates
  std::unordered_set<mesh_elem> tmp_set(std::begin(ghosted_boundary_neighbors),
  					std::end(ghosted_boundary_neighbors));
  // Convert the set to a vector
  _ghost_faces.insert(std::end(_ghost_faces),
  			   std::begin(tmp_set),std::end(tmp_set));

#ifdef USE_MPI
  SPDLOG_DEBUG("MPI Process {} has {} ghosted faces.",_comm_world.rank(), _ghost_faces.size());
#endif
}

void triangulation::shrink_local_mesh_to_owned_and_distance_neighbors()
{
  // Reset _faces to contain ONLY _local_faces and _ghost_faces.

  SPDLOG_DEBUG("Shrinking local meshes");

  // Ensure that the localized face containers have been set
  assert( _local_faces.size() != 0 );
  assert( _ghost_faces.size() != 0 );

  // remove everything from _faces and force reallocation
  _faces.clear();
  std::vector<mesh_elem>().swap(_faces);

  // Put all local and ghost faces into the _faces vector
  _faces.reserve( _local_faces.size() + _ghost_faces.size() );
  _faces.insert( _faces.end(), _local_faces.begin(), _local_faces.end() );
  _faces.insert( _faces.end(), _ghost_faces.begin(), _ghost_faces.end() );

#ifdef USE_MPI
  SPDLOG_DEBUG("MPI Process {} has {} faces after shrinking", _comm_world.rank(), _faces.size());
#endif
}



Delaunay::Vertex_handle triangulation::vertex(size_t i)
{
    return _vertexes.at(i);
}

mesh_elem triangulation::face(size_t i)
{
#if USE_MPI
    return _local_faces.at(i);
#else
    return _faces.at(i);
#endif
}

const std::vector<int>& triangulation::get_global_IDs() const
{
  return _global_IDs;
}

void triangulation::timeseries_to_file(double x, double y, std::string fname)
{
    mesh_elem m = this->find_closest_face(x, y);
    timeseries_to_file(m, fname);
}

void triangulation::timeseries_to_file(mesh_elem m, std::string fname)
{
    if (m == NULL)
    {
        CHM_THROW_EXCEPTION(mesh_error, "Couldn't find triangle at (x,y)");
    }

    m->to_file(fname);
}

void triangulation::init_vtkUnstructured_Grid(std::vector<std::string> output_variables)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    if(_write_ghost_neighbors_to_vtu)
    {
      triangles->Allocate(this->size_faces()+_ghost_faces.size());
    } else {
      triangles->Allocate(this->size_faces());
    }

    vtkSmartPointer<vtkStringArray> proj4 = vtkSmartPointer<vtkStringArray>::New();
    proj4->SetNumberOfComponents(1);
    proj4->SetName("proj4");
    proj4->InsertNextValue(_srs_wkt);

    double scale = is_geographic() == true ? 100000. : 1.;

    std::map<int, int> global_to_local_vertex_id;
    std::vector<int> global_vertex_id;

    // npoints holds the total number of points
    int npoints=0;
    for (size_t i = 0; i < this->size_faces(); i++)
    {
        mesh_elem fit = this->face(i);

        vtkSmartPointer<vtkTriangle> tri =
                vtkSmartPointer<vtkTriangle>::New();

	// loop over vertices of a face
	for (int j=0;j<3;++j){
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

    if(_write_ghost_neighbors_to_vtu)
    {

        /* Ghost neighbors */
        for (size_t i = 0; i < this->_ghost_faces.size(); i++)
        {
            mesh_elem fit = _ghost_faces[i];

            vtkSmartPointer<vtkTriangle> tri =
                    vtkSmartPointer<vtkTriangle>::New();

            // loop over vertices of a face
            for (int j=0;j<3;++j){
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

    }  // if write ghosts


    _vtk_unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    _vtk_unstructuredGrid->SetPoints(points);
    _vtk_unstructuredGrid->SetCells(VTK_TRIANGLE, triangles);
    _vtk_unstructuredGrid->GetFieldData()->AddArray(proj4);

//    vtkSmartPointer<vtkStringArray> test = vtkStringArray::SafeDownCast(_vtk_unstructuredGrid->GetFieldData()->GetAbstractArray("proj4"));
//    LOG_DEBUG << test->GetValue(0) ;


    //assume that all the faces have the same number of variables and the same types of variables
    //by this point this should be a fair assumption

    auto variables = output_variables.size() == 0 ? this->face(0)->variables() : output_variables;
    for(auto& v: variables)
    {
        data[v] = vtkSmartPointer<vtkFloatArray>::New();
        data[v]->SetName(v.c_str());
    }

    _vtu_global_id = vtkSmartPointer<vtkUnsignedLongArray>::New();
    _vtu_global_id->SetName("global_id");

    if(_write_parameters_to_vtu)
    {
        auto params = this->face(0)->parameters();
        for (auto &v: params)
        {
            data["[param] " + v] = vtkSmartPointer<vtkFloatArray>::New();
            data["[param] " + v]->SetName(("[param] " + v).c_str());
        }

        auto ics = this->face(0)->initial_conditions();
        for(auto& v: ics)
        {
            data["[ic] " + v] = vtkSmartPointer<vtkFloatArray>::New();
            data["[ic] " + v]->SetName(("[ic] " + v).c_str());
        }

        //handle elevation/aspect/slope
        data["Elevation"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Elevation"]->SetName("Elevation");

        data["Slope"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Slope"]->SetName("Slope");

        data["Aspect"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Aspect"]->SetName("Aspect");

        data["Area"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Area"]->SetName("Area");

        data["is_ghost"] = vtkSmartPointer<vtkFloatArray>::New();
        data["is_ghost"]->SetName("is_ghost");

        data["ghost_type"] = vtkSmartPointer<vtkFloatArray>::New();
        data["ghost_type"]->SetName("ghost_type");

#ifdef USE_MPI
        data["owner"] = vtkSmartPointer<vtkFloatArray>::New();
        data["owner"]->SetName("owner");
#endif
    }
    auto vec = this->face(0)->vectors();
    for(auto& v: vec)
    {
        vectors[v] = vtkSmartPointer<vtkFloatArray>::New();
        vectors[v]->SetName(v.c_str());
        vectors[v]->SetNumberOfComponents(3);
    }

    // Global vertex ids -> only need to be set here, get written in the writer
//    vertex_data["global_id"] = vtkSmartPointer<vtkFloatArray>::New();
//    vertex_data["global_id"]->SetName("global_id");
//    for(int i=0;i<npoints;++i)
//    {
//      vertex_data["global_id"]->InsertTuple1(i,global_vertex_id[i]);
//    }
}

void triangulation::init_timeseries(std::set< std::string > variables)
{
    #pragma omp parallel for
    for (size_t it = 0; it < size_faces(); it++)
    {
        auto face = this->face(it);
        face->init_time_series(variables);
    }

}

void triangulation::init_vectors(std::set<std::string>& variables)
{
#pragma omp parallel for
    for (size_t it = 0; it < size_faces(); it++)
    {
        auto face = this->face(it);
        face->init_vectors(variables);
    }
}

void triangulation::init_module_data(std::set< std::string > modules)
{
#pragma omp parallel for
    for (size_t it = 0; it < size_faces(); it++)
    {
        auto face = this->face(it);
        face->init_module_data(modules);
    }
}

void triangulation::prune_faces(std::vector<Face_handle>& faces)
{
#ifdef USE_MPI
    if(_comm_world.size() > 1)
    {
        CHM_THROW_EXCEPTION(config_error, "Cannot prune faces with ranks >1");
    }
#endif

    _faces.clear();
    _faces  = faces;

    _local_faces.clear();
    _local_faces = _faces;

    _num_faces = _num_global_faces = _faces.size(); //number of global faces

}
void triangulation::init_face_data(std::set< std::string >& timeseries,
                    std::set< std::string >& vectors,
                    std::set< std::string >& module_data)
{
    #pragma omp parallel for
        for (size_t it = 0; it < size_faces(); it++)
        {
            auto face = this->face(it);
            face->init_time_series(timeseries);
            face->init_module_data(module_data);
            face->init_vectors(vectors);
        }
	// Init data in ghost neighbors as well
	// - so they can be treated just like normal neighbors (after vars communicated)
	// - timeseries not needed here
	SPDLOG_DEBUG("######### Current _ghost_neighbors.size(): {}",_ghost_neighbors.size());
    #pragma omp parallel for
        for (size_t it = 0; it < _ghost_faces.size(); it++)
        {
            auto face = _ghost_faces.at(it);
            face->init_module_data(module_data);
            face->init_time_series(timeseries);
            face->init_vectors(vectors);
        }
}

void triangulation::update_vtk_data(std::vector<std::string> output_variables)
{
    //if we haven't inited yet, do so.
    if(!_vtk_unstructuredGrid || _terrain_deformed)
    {
        this->init_vtkUnstructured_Grid(output_variables);
    }

    auto variables = output_variables.size() == 0 ? this->face(0)->variables() : output_variables;
    auto params = this->face(0)->parameters();
    auto ics = this->face(0)->initial_conditions();
    auto vecs = this->face(0)->vectors();

    for (size_t i = 0; i < this->size_faces(); i++)
    {
        mesh_elem fit = this->face(i);

        for (auto &v: variables)
        {
            double d = (*fit)[v];
            if(d == -9999.)
            {
                d = nan("");
            }

            data[v]->InsertTuple1(i,d);
        }

        //this is mandatory now
        _vtu_global_id->InsertTuple1(i, fit->cell_global_id);

        if(_write_parameters_to_vtu)
        {
            for (auto &v: params)
            {
                double d = fit->parameter(v);
                if (d == -9999.)
                {
                    d = nan("");
                }
                data["[param] " + v]->InsertTuple1(i, d);
            }

            for (auto &v: ics)
            {
                double d = fit->get_initial_condition(v);
                if (d == -9999.)
                {
                    d = nan("");
                }
                data["[ic] " + v]->InsertTuple1(i, d);
            }

            data["Elevation"]->InsertTuple1(i,fit->get_z());
            data["Slope"]->InsertTuple1(i,fit->slope());
            data["Aspect"]->InsertTuple1(i,fit->aspect());
            data["Area"]->InsertTuple1(i,fit->get_area());
	    data["is_ghost"]->InsertTuple1(i,fit->is_ghost);
            data["ghost_type"]->InsertTuple1(i,fit->ghost_type);

#ifdef USE_MPI
	    data["owner"]->InsertTuple1(i,_comm_world.rank());
#endif
        }
        for(auto& v: vecs)
        {
            Vector_3 d = fit->face_vector(v);

            vectors[v]->InsertTuple3(i,d.x(),d.y(),d.z());
        }


    }

    if(_write_ghost_neighbors_to_vtu)
    {

    /* Ghost neighbors */
    for (size_t i = 0; i < _ghost_faces.size(); i++)
    {
        mesh_elem fit = _ghost_faces[i];

	size_t insert_offset = i + this->size_faces();


        for (auto &v: variables)
        {
            double d = -9999;
            if(fit->ghost_type == GHOST_TYPE::NEIGH)
                d = (*fit)[v];

            if(d == -9999.)
            {
                d = nan("");
            }

            data[v]->InsertTuple1(insert_offset,d);
        }

        _vtu_global_id->InsertTuple1(insert_offset,fit->cell_global_id);

        if(_write_parameters_to_vtu)
        {
            for (auto &v: params)
            {
                double d = fit->parameter(v);
                if (d == -9999.)
                {
                    d = nan("");
                }
                data["[param] " + v]->InsertTuple1(insert_offset, d);
            }

            for (auto &v: ics)
            {
                double d = fit->get_initial_condition(v);
                if (d == -9999.)
                {
                    d = nan("");
                }
                data["[ic] " + v]->InsertTuple1(insert_offset, d);
            }

            data["Elevation"]->InsertTuple1(insert_offset,fit->get_z());
            data["Slope"]->InsertTuple1(insert_offset,fit->slope());
            data["Aspect"]->InsertTuple1(insert_offset,fit->aspect());
            data["Area"]->InsertTuple1(insert_offset,fit->get_area());
	    data["is_ghost"]->InsertTuple1(insert_offset,fit->is_ghost);
            data["ghost_type"]->InsertTuple1(insert_offset,fit->ghost_type);

	    data["owner"]->InsertTuple1(insert_offset,fit->owner);
        }
        for(auto& v: vecs)
        {
            Vector_3 d = fit->face_vector(v);

            vectors[v]->InsertTuple3(insert_offset,d.x(),d.y(),d.z());
        }


    }

    } // if write ghosts

    _vtk_unstructuredGrid->GetCellData()->AddArray(_vtu_global_id);

    for(auto& m : vectors)
    {
        _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);

    }


    for(auto& m : data)
    {
        _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);
    }

    for(auto& m : vertex_data)
    {
        _vtk_unstructuredGrid->GetPointData()->AddArray(m.second);
    }

}
void triangulation::write_vtu(std::string file_name)
{
    //this now needs to be called from outside these functions
//    update_vtk_data();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file_name.c_str());
//    writer->SetCompressorType( vtkXMLUnstructuredGridWriter::CompressorType::ZLIB);
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(_vtk_unstructuredGrid);
#else
    writer->SetInputData(_vtk_unstructuredGrid);
#endif
    writer->Write();

//    write_vtp(file_name);
}

double triangulation::max_z()
{
    return _max_z;
}

double triangulation::min_z()
{
    return _min_z;
}

boost::shared_ptr<segmented_AABB> triangulation::AABB(size_t rows, size_t cols)
{
boost::shared_ptr<segmented_AABB> AABB = boost::make_shared<segmented_AABB>();
AABB->make(this,rows,cols);
return AABB;

}


void segmented_AABB::make( triangulation* domain, size_t rows, size_t cols)
{
    int n_segments = cols;
    int n_v_segments = rows;

    arma::vec midpoint_bottom_x(n_segments);
    arma::vec midpoint_top_x(n_segments);

    arma::vec midpoint_bottom_y(n_segments);
    arma::vec midpoint_top_y(n_segments);


    this->n_rows = rows;
    this->n_cols = cols;

    std::vector<K::Point_2> bboxpt;
    bboxpt.reserve( domain->number_of_vertices());

    for(auto v = domain->vertices_begin(); v!=domain->vertices_end();v++ )
    {
        bboxpt.push_back(K::Point_2(v->point().x(),v->point().y()));
    }

    K::Iso_rectangle_2 rot_bbox = CGAL::bounding_box(bboxpt.begin(), bboxpt.end());

    //need to construct new axis aligned BBR
//    arma::u32 index;

    //left most pt
    double x_left = rot_bbox.vertex(0).x();
//    double y_left = rot_bbox.vertex(0).y();

    //right most pt
    double x_right = rot_bbox.vertex(2).x();
//    double y_right = rot_bbox.vertex(2).y();

    //bottom most pt
    double y_bottom = rot_bbox.vertex(1).y();
//    double x_bottom = rot_bbox.vertex(1).x();

    //top most pt
    double y_top =rot_bbox.vertex(3).y();
//    double x_top = rot_bbox.vertex(3).x();


    //horizontal step size
    double h_dx = (x_right - x_left) / n_segments;
    double v_dy = (y_top - y_bottom) / n_v_segments;


    //resize the row dimension
    m_grid.resize(n_v_segments);

    //start top left
    //loop over rows
    for (int i = 0; i < n_v_segments; i++)
    {
        //set the y coordinate to the current row top left coordinate
        double h_y = y_top - v_dy*i;
        double h_x = x_left;

        //set number of cols
        m_grid[i].resize(n_segments);


        //loop over columns
        for (int j = 0; j < n_segments; j++)
        {

            arma::mat* t = new arma::mat(5, 2);


            *t << h_x << h_y - v_dy << arma::endr // bottom left
                    << h_x + h_dx << h_y - v_dy << arma::endr //bottom right
                    << h_x + h_dx << h_y << arma::endr // top right
                    << h_x << h_y << arma::endr //top left
                    << h_x << h_y - v_dy << arma::endr; // bottom left


            m_grid[i][j] = new rect(t);
//            m_grid[i][j]->triangles.reserve( ceil(domain->number_of_vertices()/(rows*cols)) ); //estimate an approx number of triangles/rect. Try to cut down on repeat memory alloc.
            h_x = h_x + h_dx;
        }
    }


}

rect* segmented_AABB::get_rect(size_t row, size_t col)
{
    return m_grid[row][col];
}

segmented_AABB::segmented_AABB()
{
    n_rows = 0;
    n_cols = 0;
}

bool segmented_AABB::pt_in_rect(Delaunay::Vertex_handle v, rect* r)
{
    return pt_in_rect(v->point().x(), v->point().y(), r);
}

bool segmented_AABB::pt_in_rect(double x, double y, rect* r)
{
    for (int i = 0; i < 4; i++)
    {
        double x0 = r->coord->operator()(i, 0);
        double y0 = r->coord->operator()(i, 1);
        double x1 = r->coord->operator()(i + 1, 0);
        double y1 = r->coord->operator()(i + 1, 1);

        double pt = ((y - y0)*(x1 - x0) - (x - x0)*(y1 - y0));
        if (pt < 0)
            return false;
        else if (pt == 0)
            return false;
    }
    return true;
}

segmented_AABB::~segmented_AABB()
{
    for (auto it = m_grid.begin(); it != m_grid.end(); it++)
    {
        for (auto jt = it->begin(); jt != it->end(); jt++)
        {
            delete *jt;
        }
    }

}
	// compute_balanced_partition_sizes();
