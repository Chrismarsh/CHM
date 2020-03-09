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
 //   LOG_WARNING << "No Matlab engine, plotting and all Matlab functionality will be disabled";
#ifdef MATLAB
    _engine = NULL;
    _gfx = NULL;
#endif
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

#endif
}

#ifdef MATLAB

triangulation::triangulation(boost::shared_ptr<maw::matlab_engine> engine)
{
    _engine = engine;
    _gfx = boost::make_shared<maw::graphics>(_engine.get());
    _size = 0;
}
#endif

triangulation::~triangulation()
{
#ifdef MATLAB
    _engine = NULL;
    _gfx = NULL;
#endif
}

//void triangulation::init(vector x, vector y, vector z)
//{
//    //    _mesh = boost::make_shared<Delaunay>();
//    for (size_t i = 0; i < x->size(); i++)
//    {
//        Point_3 p(x->at(i), y->at(i), z->at(i));
//        Delaunay::Vertex_handle Vh = this->insert(p);
//        Vh->set_id(i);
//        ++i;
//    }
//
//    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
//            fit != this->finite_faces_end(); ++fit)
//    {
//        mesh_elem face = fit;
//        _faces.push_back(face);
//    }
//
//    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(this->size_faces()) + " triangles";
//}

std::string triangulation::proj4()
{
    return _srs_wkt;
}

void triangulation::write_param_to_vtu(bool write_param)
{
    _write_parameters_to_vtu = write_param;
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
      if(!f->_is_ghost &&
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

#ifdef NOMATLAB

void triangulation::plot_time_series(double x, double y, std::string ID)
{
    if (!_engine)
    {
        LOG_WARNING << "No Matlab engine, plotting is disabled";
        return;
    }
    mesh_elem m = this->locate_face(x, y);

    if (m == NULL)
        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Couldn't find triangle at (x,y)"));

    maw::d_vec v(new arma::vec(m->face_time_series(ID).size()));

    timeseries::variable_vec ts = m->face_time_series(ID);

    for (size_t i = 0; i < v->size(); i++)
    {
        (*v)(i) = ts.at(i);
    }

    _engine->put_double_matrix(ID, v);
    double handle = _gfx->plot_line(ID);
    _gfx->add_title(ID);
    _gfx->spin_until_close(handle);
}
#endif

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



    size_t nvertex_toread = mesh.get<size_t>("mesh.nvertex");
    LOG_DEBUG << "Reading in #vertex=" << nvertex_toread;
    size_t i=0;

    //paraview struggles with lat/long as it doesn't seem to have the accuracy. So we need to scale up lat-long.
    int is_geographic = mesh.get<int>("mesh.is_geographic");
    if( is_geographic == 1)
    {
        _is_geographic = true;
    }
    _srs_wkt = mesh.get<std::string>("mesh.proj4","");


    if(_srs_wkt == "")
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info("proj4 field in .mesh file is empty!"));
    }

    for (auto &itr : mesh.get_child("mesh.vertex"))
    {
        std::vector<double> items;
        //iterate over the vertex triples
        for(auto& jtr: itr.second)
        {
            items.push_back(jtr.second.get_value<double>());
        }
        Point_3 pt( items[0], items[1], items[2]);

        _max_z = std::max(_max_z,items[2]);
        _min_z = std::min(_min_z,items[2]);

        Vertex_handle Vh = this->create_vertex();
        Vh->set_point(pt);
        Vh->set_id(i);
        _vertexes.push_back(Vh);
        i++;
    }
    _num_vertex = this->number_of_vertices();
    LOG_DEBUG << "# nodes created = " << _num_vertex;

    if( this->number_of_vertices() != nvertex_toread)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                "Expected: " + std::to_string(nvertex_toread) + " vertex, got: " + std::to_string(this->number_of_vertices())));
    }

    //read in faces
    size_t num_elem = mesh.get<int>("mesh.nelem");
    LOG_DEBUG << "Reading in #elem = " << num_elem;
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

// If we aren't using MPI, we can build the centre points here for efficiency.
// these centre points will be used to build the spatial search tree. If we are using MPI,
// we need to wait until we've figured out the per-node triangle partition so we can build
// a per-node spatial search tree that only takes into account this node's elements.
// If MPI, we do that at the end of this function

#ifndef USE_MPI
        Point_2 pt2(face->center().x(),face->center().y());
        center_points.push_back(pt2);
#endif

    }

#ifndef USE_MPI
    //make the search tree
    dD_tree = boost::make_shared<Tree>(boost::make_zip_iterator(boost::make_tuple( center_points.begin(),_faces.begin() )),
                                       boost::make_zip_iterator(boost::make_tuple( center_points.end(), _faces.end() ) )
    );
#endif

    _num_faces = this->number_of_faces();

    LOG_DEBUG << "Created a mesh with " << this->size_faces() << " triangles";

    if( this->number_of_faces() != num_elem)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                "Expected: " + std::to_string(num_elem) + " elems, got: " + std::to_string(this->size_faces())));
    }

    LOG_DEBUG << "Building face neighbours";
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
            BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                    "Face " + std::to_string(i) + " has out of bound neighbours."));
        }
        //-1 is now the no neighbour value
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
                LOG_WARNING << "Parameter " + name + " is zero length and will be ignored.";
            } else {
                _parameters.insert(name);
            }

        }

        // init the storage, which builds the mphf
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

            LOG_DEBUG << "Applying parameter: " << name;

            for (auto &jtr : itr.second)
            {
                auto face = _faces.at(i);
                auto value = jtr.second.get_value<double>();
                face->parameter(name) = value;
                i++;
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
            LOG_DEBUG << "Applying IC: " << name;
            for (auto &jtr : itr.second)
            {
                auto face = _faces.at(i);
                double value = jtr.second.get_value<double>();
    //            alue == -9999. ? value = nan("") : value;
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
        LOG_DEBUG << "No face permutation.";
    }

    partition_mesh();
#ifdef USE_MPI
    _num_faces = _local_faces.size();
    size_t total_num_faces = _faces.size();
    determine_local_boundary_faces();
    determine_process_ghost_faces_nearest_neighbours();

    // should make this parallel
    for(size_t ii=0; ii < total_num_faces; ++ii)
    {
        auto face = _faces.at(ii);
        Point_2 pt2(face->center().x(),face->center().y());
        center_points.push_back(pt2);

    }
    //make the search tree
    dD_tree = boost::make_shared<Tree>(boost::make_zip_iterator(boost::make_tuple( center_points.begin(),_faces.begin() )),
                                       boost::make_zip_iterator(boost::make_tuple( center_points.end(),  _faces.end() ) )
				       );

#endif // USE_MPI

    // determining ghost faces requires the dD_tree to be set up
    // determine_process_ghost_faces_by_distance(100.);

    // shrink the local mesh
    // shrink_local_mesh_to_owned_and_distance_neighbours();

  std::vector<double> temp_slope(_num_faces);

#pragma omp parallel for
  for (size_t i = 0; i < _num_faces; i++)
  {

    auto f = face(i);
    std::vector<boost::tuple<double, double, double> > u;
    for (size_t j = 0; j < 3; j++)
    {
      auto neigh = f->neighbor(j);
      if (neigh != nullptr && !neigh->_is_ghost)
        u.push_back(boost::make_tuple(neigh->get_x(), neigh->get_y(), neigh->slope()));
    }

    auto query = boost::make_tuple(f->get_x(), f->get_y(), f->get_z());

    interpolation interp(interp_alg::tpspline);
    double new_slope = f->slope();

    if(u.size() > 0)
    {
      new_slope = interp(u, query);
    }

    temp_slope.at(i) = new_slope;
  }

#pragma omp parallel for
  for (size_t i = 0; i < size_faces(); i++)
  {
    auto f = face(i);
    f->_slope = temp_slope.at(i);
    //init these
    f->aspect();
    f->center();
    f->normal();
  }

    // TODO need to re-setup dD_tree to only consider the _faces after shrinking
    // -  Note this likely allows us to remove the the ifndef USE_MPI from earlier in this routine,
    //    as we have to do a global pass and a local pass anyway

}

void triangulation::reorder_faces(std::vector<size_t> permutation)
{
  // NOTE: be careful with evaluating this, the 'cell_global_id's and a
  // cell's position in the '_faces' vec are unrelated. They must both
  // be modified (ie. renumber the 'cell_global_id's, AND sort the '_faces'
  // vector) before the new ordering is consistent.

  assert( permutation.size() == size_faces() );
  LOG_DEBUG << "Reordering faces";



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

void triangulation::partition_mesh()
{
  /*
    Basically, this first attempt splits the number of mesh points as equally as possible
  */
  // 1. Determine number of (processor) locally owned faces
  // 2. Determine (processor) locally owned indices of faces
  // 3. Set locally owned _is_ghost=false



  size_t total_num_faces = _faces.size();

#ifdef USE_MPI

  LOG_DEBUG << "Partitioning mesh";

  // Set up so that all processors know how 'big' all other processors are
  std::vector<int> num_faces_in_partition(_comm_world.size(),
  					  total_num_faces/_comm_world.size());
  for (unsigned int i=0;i<total_num_faces%_comm_world.size();++i) {
    num_faces_in_partition[i]++;
  }

  // each processor only knows its own start and end indices
  size_t face_start_idx = 0;
  size_t face_end_idx = num_faces_in_partition[0]-1;
  for (int i=1;i<=_comm_world.rank();++i) {
    face_start_idx += num_faces_in_partition[i-1];
    face_end_idx += num_faces_in_partition[i];
  }


  // Set size of vector containing locally owned faces
  _local_faces.resize(num_faces_in_partition[_comm_world.rank()]);

#pragma omp parallel for
  for(int local_ind=0;local_ind<_local_faces.size();++local_ind)
  {

	     size_t global_ind = face_start_idx + local_ind;
	     _faces.at(global_ind)->_is_ghost = false;
	     _faces.at(global_ind)->cell_local_id = local_ind;
	     _local_faces[local_ind] = _faces.at(global_ind);

  }


  LOG_DEBUG << "MPI Process " << _comm_world.rank() << ": start " << face_start_idx << ", end " << face_end_idx << ", number " << _local_faces.size();


#else // do not USE_MPI

#pragma omp parallel for
  for(size_t i=0;i<total_num_faces;++i)
  {
    _faces.at(i)->_is_ghost = false;
    _faces.at(i)->cell_local_id = i; // Mesh has been (potentially) reordered before this point. Set the local_id correctly
  }
  LOG_DEBUG << "Face numbering : start 0, end " << (total_num_faces-1) << ", number " << _local_faces.size();

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

		 int num_owned_neighbours = 0;
		 for (int neigh_index = 0; neigh_index < 3; ++neigh_index)
		   {

		     auto neigh = face->neighbor(neigh_index);

		     // Test status of neighbour
		     if (neigh == nullptr)
		       {
			 th_local_boundary_faces[omp_get_thread_num()].push_back(std::make_pair(face,true));
			 num_owned_neighbours=3; // set this to avoid triggering the post-loop if statement
			 break;
		       } else
		       {
			 if (neigh->_is_ghost == false)
			   {
			     num_owned_neighbours++;
			   }
		       }
		   }

		 // If we don't own 3 neighbours, we are a local, but not a global boundary face
		 if( num_owned_neighbours<3 ) {
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

#endif
}

void triangulation::determine_process_ghost_faces_nearest_neighbours()
{
  // NOTE that this algorithm is not implemented for multithread
  // - multithread can be implemented similarly to determine_local_boundary_faces
  //   - not a priority since number of boundary faces should be small

  // Ensure that the local boundary faces have been determined, but ghost nearest neighbours have not been set
  assert( _boundary_faces.size() != 0 );
  assert( _ghost_neighbours.size() == 0 );

  LOG_DEBUG << "Determining ghost region info";

  // Vector for append speed
  std::vector< mesh_elem > ghosted_boundary_nearest_neighbours;

  for(size_t face_index=0; face_index< _boundary_faces.size(); ++face_index)
  {
    // face_index is a local index... get the face handle
    auto face = _boundary_faces.at(face_index).first;
    // append the ghosted nearest neighbours
    for(int i = 0; i < 3; ++i)
    {
        auto neigh = face->neighbor(i);
        if(neigh != nullptr && neigh->_is_ghost)
            ghosted_boundary_nearest_neighbours.push_back(neigh);
    }
  }

  // Convert to a set to remove duplicates
  std::unordered_set<mesh_elem> tmp_set(std::begin(ghosted_boundary_nearest_neighbours),
					std::end(ghosted_boundary_nearest_neighbours));
  // Convert the set to a vector
  _ghost_neighbours.insert(std::end(_ghost_faces),
			   std::begin(tmp_set),std::end(tmp_set));
#ifdef USE_MPI
  LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _ghost_neighbours.size() << " ghosted nearest neighbours.";
#endif
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

  // otherwise, visit face, and move on to neighbours
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

  // Ensure that the local boundary faces have been determined, but ghost neighbours have not been set
  assert( _boundary_faces.size() != 0 );
  assert( _ghost_faces.size() == 0 );

  // Vector for append speed
  std::vector< mesh_elem > ghosted_boundary_neighbours;

  for(size_t face_index=0; face_index< _boundary_faces.size(); ++face_index)
  {
    // face_index is a local index... get the face handle
    auto face = _boundary_faces.at(face_index).first;
    // separate the ghosted and non-ghosted neighbours
    // std::vector<mesh_elem> current_neighbours = find_faces_in_radius(face->center().x(),face->center().y(), max_distance);
    std::vector<mesh_elem> current_neighbours = dfs_to_max_distance(face, max_distance);
    auto pivot = std::partition(std::begin(current_neighbours),std::end(current_neighbours),
    				[] (mesh_elem neigh) {
    				  return neigh->_is_ghost == true;
    				});
    current_neighbours.erase(pivot,std::end(current_neighbours));

    ghosted_boundary_neighbours.insert(std::end(ghosted_boundary_neighbours),
    				       current_neighbours.begin(),current_neighbours.end());
  }

  // Convert to a set to remove duplicates
  std::unordered_set<mesh_elem> tmp_set(std::begin(ghosted_boundary_neighbours),
  					std::end(ghosted_boundary_neighbours));
  // Convert the set to a vector
  _ghost_faces.insert(std::end(_ghost_faces),
  			   std::begin(tmp_set),std::end(tmp_set));

#ifdef USE_MPI
  LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _ghost_faces.size() << " ghosted faces.";
#endif
}

void triangulation::shrink_local_mesh_to_owned_and_distance_neighbours()
{
  // Reset _faces to contain ONLY _local_faces and _ghost_faces.

  LOG_DEBUG << "Shrinking local meshes";

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
  LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _faces.size() << " faces after shrinking";
#endif
}

void triangulation::populate_face_station_lists()
{

  LOG_DEBUG << "Populating each face's station list";

  for (auto& f : _faces)
    {
      if ( f->_stations.size() == 0 ){
	auto stations = _global->get_stations(f->get_x(), f->get_y());
	f->_stations.insert(std::end(f->_stations), std::begin(stations), std::end(stations));
      }  else {
	BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Face station list already populated."));
      }
    }

}

void triangulation::populate_distributed_station_lists()
{

  assert( _stations.size() == 0 ); // A mesh's stationlist must not have already been set

  using th_safe_multicontainer_type = std::vector< std::shared_ptr<station> >[];
  std::unique_ptr< th_safe_multicontainer_type > th_local_stations;
  std::vector< std::shared_ptr<station> > mpi_local_stations;

  LOG_DEBUG << "Populating each MPI process's station list";

#pragma omp parallel
  {
    // We want an array of vectors, so that OMP threads can increment them
    // separately, then join them afterwards
#pragma omp single
    {
      th_local_stations = std::make_unique< th_safe_multicontainer_type >(omp_get_num_threads());
    }
#pragma omp for
    for(size_t face_index=0; face_index< _local_faces.size(); ++face_index)
      {
	// face_index is a local index... get the face handle
	auto face = _local_faces.at(face_index);
	if ( face->stations().size() == 0 ){ // only perform if faces' stationlists are set
	  BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Face station lists must be populated before populating distributed MPI station lists."));
	}
	for (auto &p : face->stations())
	  {
	    th_local_stations[omp_get_thread_num()].push_back(p);
	  }
      }
    // Join the vectors via a single thread in t operations
    //  NOTE future optimizations:
    //   - reserve space for insertions into mpi_local_stations
    //   - can be done recursively in log2(num_threads) operations
#pragma omp single
    {
      for(int thread_idx=0;thread_idx<omp_get_num_threads();++thread_idx)
	{
	  mpi_local_stations.insert(std::end(mpi_local_stations),
				    std::begin(th_local_stations[thread_idx]), std::end(th_local_stations[thread_idx]));
	}
    }

  }

  // Remove duplicates by converting to a set
  std::unordered_set< std::shared_ptr<station> > tmp_set(std::begin(mpi_local_stations), std::end(mpi_local_stations));

  // Store the local stations in the triangulations mpi-local stationslist vector
  _stations.insert( std::end(_stations), std::begin(tmp_set), std::end(tmp_set));
#ifdef USE_MPI
  LOG_DEBUG << "MPI Process " << _comm_world.rank() << " has " << _stations.size() << " locally owned stations.";
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

void triangulation::timeseries_to_file(double x, double y, std::string fname)
{
    mesh_elem m = this->find_closest_face(x, y);
    timeseries_to_file(m, fname);
}

void triangulation::timeseries_to_file(mesh_elem m, std::string fname)
{
    if (m == NULL)
        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Couldn't find triangle at (x,y)"));

    m->to_file(fname);
}
#ifdef NOMATLAB

void triangulation::plot(std::string ID)
{
    if (!_engine)
    {
        LOG_WARNING << "No Matlab engine, plotting is disabled";
        return;
    }
    LOG_DEBUG << "Sending triangulation to matlab...";



    maw::d_mat tri(new arma::mat(_size, 3));
    maw::d_mat xyz(new arma::mat(_data_size, 3));
    maw::d_mat cdata(new arma::mat(_size, 1));

    size_t i = 0;

    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
            fit != this->finite_faces_end(); ++fit)
    {
        auto face = fit;

        (*tri)(i, 0) = face->vertex(0)->get_id() + 1; //+1 because matlab indexing starts at 1
        (*tri)(i, 1) = face->vertex(1)->get_id() + 1;
        (*tri)(i, 2) = face->vertex(2)->get_id() + 1;

        Delaunay::Triangle t = this->triangle(face);

        //        std::cout  << "xyz1:" << (*tri)(i, 0)-1 <<std::endl;
        (*xyz)((*tri)(i, 0) - 1, 0) = t[0].x();
        (*xyz)((*tri)(i, 0) - 1, 1) = t[0].y();
        (*xyz)((*tri)(i, 0) - 1, 2) = t[0].z();


        //        std::cout  << "xyz2:"<<(*tri)(i, 1)-1 <<std::endl;
        (*xyz)((*tri)(i, 1) - 1, 0) = t[1].x();
        (*xyz)((*tri)(i, 1) - 1, 1) = t[1].y();
        (*xyz)((*tri)(i, 1) - 1, 2) = t[1].z();

        //        std::cout  << "xyz3:"<<(*tri)(i, 2)-1 <<std::endl;
        (*xyz)((*tri)(i, 2) - 1, 0) = t[2].x();
        (*xyz)((*tri)(i, 2) - 1, 1) = t[2].y();
        (*xyz)((*tri)(i, 2) - 1, 2) = t[2].z();


        double d = fit->face_data(ID);
        (*cdata)(i) = d;
        //        std::cout << i <<std::endl;
        ++i;
    }


    LOG_DEBUG << "Sending data to matlab...";
    _engine->put_double_matrix("tri", tri);
    _engine->put_double_matrix("elevation_data", xyz);
    _engine->put_double_matrix("face_data", cdata);

    double handle = _gfx->plot_patch("[elevation_data(:,1) elevation_data(:,2) elevation_data(:,3)]", "tri", "face_data(:)");
    _gfx->add_title(ID);

    _gfx->spin_until_close(handle);
    //    _engine->evaluate("save lol.mat");
    _engine->evaluate("clear tri elevation_data face_data");
}
#endif

void triangulation::init_vtkUnstructured_Grid(std::vector<std::string> output_variables)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(this->_num_vertex);

    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->Allocate(this->_num_vertex);

    vtkSmartPointer<vtkStringArray> proj4 = vtkSmartPointer<vtkStringArray>::New();
    proj4->SetNumberOfComponents(1);
    proj4->SetName("proj4");
    proj4->InsertNextValue(_srs_wkt);

    double scale = is_geographic() == true ? 100000. : 1.;

    for (size_t i = 0; i < this->size_faces(); i++)
    {
        mesh_elem fit = this->face(i);

        vtkSmartPointer<vtkTriangle> tri =
                vtkSmartPointer<vtkTriangle>::New();

        tri->GetPointIds()->SetId(0, fit->vertex(0)->get_id());
        tri->GetPointIds()->SetId(1, fit->vertex(1)->get_id());
        tri->GetPointIds()->SetId(2, fit->vertex(2)->get_id());

        points->SetPoint(fit->vertex(0)->get_id(), face(i)->vertex(0)->point().x()*scale, face(i)->vertex(0)->point().y()*scale, face(i)->vertex(0)->point().z());
        points->SetPoint(fit->vertex(1)->get_id(), face(i)->vertex(1)->point().x()*scale, face(i)->vertex(1)->point().y()*scale, face(i)->vertex(1)->point().z());
        points->SetPoint(fit->vertex(2)->get_id(), face(i)->vertex(2)->point().x()*scale, face(i)->vertex(2)->point().y()*scale, face(i)->vertex(2)->point().z());

        triangles->InsertNextCell(tri);
    }

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
    }
    auto vec = this->face(0)->vectors();
    for(auto& v: vec)
    {
        vectors[v] = vtkSmartPointer<vtkFloatArray>::New();
        vectors[v]->SetName(v.c_str());
        vectors[v]->SetNumberOfComponents(3);
    }
}

void triangulation::init_timeseries(std::set< std::string > variables)
{
    #pragma omp parallel for
    for (size_t it = 0; it < _mesh->size_faces(); it++)
    {
        auto face = _mesh->face(it);
        face->init_time_series(variables);
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
        }
        for(auto& v: vecs)
        {
            Vector_3 d = fit->face_vector(v);

            vectors[v]->InsertTuple3(i,d.x(),d.y(),d.z());
        }


    }

    for(auto& m : vectors)
    {
        _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);

    }


    for(auto& m : data)
    {
        _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);
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
