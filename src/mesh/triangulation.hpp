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

#include "interpolation.hpp"


#include "station.hpp"
#include "global.hpp"

//for valgrind, remove
#define CGAL_DISABLE_ROUNDING_MATH_CHECK
// spatial tree  -- cgal includes
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Splitters.h>
#include <CGAL/Euclidean_distance.h>

#ifdef _OPENMP
#include <omp.h>
#else
// we need to define these and have them return constant values under a no-omp situation
inline int omp_get_num_threads() { return 1;}
inline int omp_get_thread_num() { return 0;}
inline int omp_get_max_threads() { return 1;}
#endif

#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <vector>
#include <set>
#include <unordered_set>
#include <stack>
#include <fstream>
#include <utility>


#include <armadillo>

#include <ogr_spatialref.h>

#ifdef USE_SPARSEHASH
#include <sparsehash/dense_hash_map>
#else
#include <unordered_map>
#endif

// hash functions
#include "utility/BBhash.h"
#include "utility/wyhash.h"

// json reader
#include "utility/readjson.hpp"

// boost includes
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_ds_face_base_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Triangulation_2.h>


#include <tbb/concurrent_vector.h>
#include <tbb/parallel_sort.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

namespace pt = boost::property_tree;

//required for the spatial searching
#include <CGAL/Search_traits_2.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>

#include <CGAL/Euclidean_distance.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>


#ifdef NOMATLAB
#include "libmaw.h"
#endif

#ifdef USE_MPI
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#endif

#include "vertex.hpp"
#include "timeseries.hpp"
#include "math/coordinates.hpp"
#include "utility/xxh64.hpp"

#include "timeseries/variablestorage.hpp"

// #include "hdf5.h"
#include "H5Cpp.h"
using namespace H5;
/**
* \struct face_info
* A way of embedding arbirtrary data into the face. This is how modules should store their data.
*/
struct face_info
{

    virtual ~face_info()
    {
    };
};

//fwd decl
class segmented_AABB;
class triangulation;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Triangle_3 Triangle_3;
typedef K::Point_3 Point_3;
typedef K::Point_2 Point_2;
typedef K::Vector_3 Vector_3;
typedef K::Vector_2 Vector_2;
typedef CGAL::Projection_traits_xy_3<K> Gt; //allows for using 2D algorithms on the 3D points

typedef ex_vertex<Gt> Vb; //custom vertex class




/**
* \class face
* \brief Defines the triangle face
*/

template < class Gt, class Fb = CGAL::Constrained_triangulation_face_base_2<Gt> >
class face : public Fb
{
    friend class triangulation;
public:

    /**
     * Multiuse flag for algorithms to use.
     * Does not make any gaurantee about state.
     */
    bool coloured;

//    face_info* info;

    /**
    * \typedef Vertex_handle Handle to a vertex
    */
    typedef typename Fb::Vertex_handle Vertex_handle;

    /**
    * \typedef Face_handle Handle to a face
    */
    typedef typename Fb::Face_handle Face_handle;


    template < typename TDS2 >
    struct Rebind_TDS
    {
        typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
        typedef face<Gt, Fb2> Other;
    };

    face();
    face(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2);
    face(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2);
    face(Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Face_handle n0, Face_handle n1, Face_handle n2, bool c0, bool c1, bool c2);

    ~face();

    /**
    * Aspect of the face. North = 0, CW . Calculated on first use, subsequent usages will not recalculate
    * \return Face aspect [rad]
    */
    double aspect();

    /**
    * Slope of the face. Calculated on first use, subsequent usages will not recalculate
    * \return slope [rad]
    */
    double slope();

    /**
    * Normalized face normal. Calculated on first use, subsequent usages will not recalculate
    */
    Vector_3 normal();

    /**
    * Center of the face as defined by a centroid. Calculated on first use, subsequent usages will not recalculate
    */
    Point_3 center();

    /**
     * Exactly the same as find_closest_face in triangulation but uses the current face's center
     * @param azimuth
     * @param distance
     * @return
     */
    const Face_handle find_closest_face(double azimuth, double distance);

    /**
     * Returns the ith edge's length. Refering to the docs here
     * http://doc.cgal.org/latest/Triangulation_2/classCGAL_1_1Triangulation__2.html
     * edgth_length(i) returns the length of the edge shared with neighbor(i),
     * opposite vertex(i)
     * @param i
     * @return
     */
    double edge_length(int i);

    /**
     * Returns a 2D vector that defines the edge i.
     * The Edge(i) is edge common to faces f and f.neighbor(i). It is also the edge joining the vertices vertex(cw(i)) and vertex(ccw(i)) of f.
     * 0-3, ccw
     * @param i
     * @return
     */
    Vector_2 edge(int i);

    /**
     * Returns the midpoint (x,y) of edge i
     * @param i
     * @return
     */
    Point_2 edge_midpoint(int i);

    /**
     * Returns the verticies a and b that make up edge i
     * @param i
     * @return
     */
    std::pair<Point_2,Point_2> edge_vertexes(int i);

    /**
     * Returns the 2D unit normal to a triangle edge i faceing outward.
     * Edge i corresponds to the same edge as edge_length.
     * @param i
     * @return
     */
    Vector_2 edge_unit_normal(int i);

    /// Returns the nearest station to the face
    /// @return
    std::shared_ptr<station>& nearest_station();

    /**
    * Returns the face's vector of stations
    */
  std::vector<std::shared_ptr<station>>& stations();

    /**
    * Checks if a point x,y is within the face
    * \return true if this face contains the point x,y
    */
    bool contains(double x, double y);

    /**
    * Checks if a point p is within the face
    * \return true if this face contains the point p
    */
    bool contains(Point_3 p);

    /**
    * Checks if the given face intersects this face in 2D. The algorithm checks if the specified face contains any of the 3 vertexs or center of this face.
    * \param fh Face handle to the other face
    * \return true if this face and fh intersects
    */
    bool intersects(Face_handle fh);

    bool has_vegetation();

    double veg_attribute(const std::string &variable);

    /**
     * Sets the vector for the given variable.
     * Does not support timeseries output.
     * There are no checks on inter-module vector guarantees like normal variables
     * Proceed with caution
     * @param variable
     * @param v
     */
    void set_face_vector(const std::string& variable, Vector_3 v);

    double& operator[](const uint64_t& variable);
    double& operator[](const std::string& variable);
    /**
     * Returns the face vector for a specified variable
     * @param variable
     * @return
     */
    Vector_3 face_vector(const std::string& variable);

    /**
    * Initializes  this faces variable storage
    * \param variables Names of the variables to add
    */
    void init_time_series(std::set<std::string>& variables);

    /**
    * Initializes  this faces vector storage
    * \param variables Names of the vectors to add
    */
    void init_vectors(std::set<std::string>& variables);

    /**
        * Initializes  this face's parameter storage
        * \param variables Names of the variables to add
    */
    void init_parameters(std::set<std::string>& parameters);

    /**
    * Initializes this face's module data storage
    * \param variables Names of the modules that will store data
*/
    void init_module_data(std::set<std::string>& modules);

    /**
    * Obtains the timeseries associated with the given variable
    * \param ID variable
    */
    timeseries::variable_vec face_time_series(std::string ID);

    /**
     * Returns true if the variable is present
     */
    bool has(const std::string& variable);

    /**
    * Returns true if the variable is present
    */
    bool has(const uint64_t& hash);

    /**
     * Returns a list of variables in this face's timeseries
     * \return Vector of variable names
     */
    std::vector<std::string> variables();

    /**
    * Returns the underlying timeseries object
    * \return Pointer to the underlying timeseries
    */
    boost::shared_ptr<timeseries> get_underlying_timeseries();

    /**
    * Returns the iterator of the current timestep.
    */
    timeseries::iterator now();

    /**
     * Increments the face to the next timestep
     */
    void next();

    /**
     * Resets the face to the start of the timeseries
     */
    void reset_to_begining();

    /**
     * Returns the xcoordinates of the center of the face
     * @return x UTM coordinate
     */
    double get_x();

    /**
     * Returns the ycoordinates of the center of the face
     * @return y UTM coordinate
     */
    double get_y();

    /**
     * Returns the zcoordinates of the center of the face
     * @return  elevation
     */
    double get_z();

    /**
     * Returns the sub-triangle z coordinates of the Point 2
     * @return  elevation
     */
    double get_subgrid_z(Point_2 query);


    /**
     * Get triangle area (m)
     */
    double get_area();
    /**
     * Saves this face's timeseries to a file
     * @param fname specified file name
     */
    void to_file(std::string fname);

    template<typename T>
    T& get_module_data(const std::string &module);

//    void set_module_data(const std::string &module, face_info *fi);

    template<typename T>
    T& make_module_data(const std::string &module);

    std::string _debug_name; //for debugging to find the elem that we want
    int _debug_ID; //also for debugging. ID == the position in the output order, starting at 0
    size_t cell_global_id;
    size_t cell_local_id;


    /**
     * Gets the face parameter value. E.g., landcover type
     * @param key
     * @return
     */
    double& parameter(const uint64_t& hash);
    double& parameter(const std::string& variable);
    /**
     * Sets the parameter on the face to the given value. Parameter doesn't have to exist. Do not use to store model output.
     * @param key
     * @param value
     */
//    void set_parameter(std::string key, double value);


    std::vector<std::string> parameters();
    bool has_parameter(const uint64_t& hash);
    bool has_parameter(const std::string& variable);

    void set_initial_condition(std::string key,double value);
    double get_initial_condition(std::string key);

    /**
     * Returns a vector of names of all intial conditions
     * @return
     */
    std::vector<std::string> initial_conditions();

    /**
     * Returns a vector of all the names of all the vector data
     * @return
     */
    std::vector<std::string> vectors() ;

    bool has_initial_condition(std::string key);

    bool _is_geographic;

    bool is_ghost=true;
    int ghost_type;

    int  owner;  // MPI process that owns the face

private:

    OGRSpatialReference _insrs; // will hold the crs of the mesh
    OGRSpatialReference _face_utm_srs; // will hold the crs of the face,


    double _slope;
    double _azimuth;
    double _area;
    double _x;
    double _y;
    double _z;


    //hold a pointer *back* to the triangulation. This let's use query triangles at distance X, etc
    //that allows for using data::parallel modules w/o having to use domain parallel.
    //const so we can't modify the domain via this as thar be dragons
    triangulation* _domain;

    boost::shared_ptr<Point_3> _center;
    boost::shared_ptr<Vector_3> _normal;


    variablestorage<double> _variables;
    variablestorage<double> _parameters;

    variablestorage< std::unique_ptr<face_info>> _module_face_data;
    variablestorage<double> _initial_conditions;
    variablestorage< Vector_3> _module_face_vectors; //holds vector components, currently no checks on anything. Proceed with caution.

    boost::shared_ptr<timeseries> _data;
    timeseries::iterator _itr;

    std::vector<std::shared_ptr<station>> _stations;

    std::shared_ptr<station> _nearest_station;

};

typedef face<Gt> Fb; //custom face class
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Delaunay;

typedef Delaunay::Face_handle mesh_elem;
typedef boost::shared_ptr<tbb::concurrent_vector<double>  > vector;

//search tree typedefs
//http://doc.cgal.org/latest/Spatial_searching/index.html
typedef K::Point_2 Point_2;
typedef boost::tuple<Point_2, mesh_elem > Point_and_face;
typedef CGAL::Search_traits_2<K>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_face,
            CGAL::Nth_of_tuple_property_map<0, Point_and_face>,
            Traits_base>                                              Traits;

//Sliding_midpoint spatial search tree. Better stability for searching with the coordinate systems we use
typedef CGAL::Sliding_midpoint<Traits> Splitter;
typedef CGAL::Kd_tree<Traits,Splitter> Tree;

// for nearest face search
typedef CGAL::Orthogonal_k_neighbor_search <
        Traits,
        typename CGAL::internal::Spatial_searching_default_distance<Traits>::type,
        Splitter > K_neighbor_search;

// used for faces in radius
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;


/**
*
*/
class triangulation
: public Delaunay
{
public:
    triangulation();
#ifdef MATLAB
/**
* If the Matlab engine integration is enabled, a pointer to the engine instance must be passed in.
* \param engine Pointer to a valid matlab engine instant
*/
    triangulation(boost::shared_ptr<maw::matlab_engine> engine);
#endif
    ~triangulation();



    /**
    * Creates an axis alligned bounding box that is segmented into rows x cols
    * \param rows number of rows in the AABB
    * \param cols number of cols in the AABB
    */
    boost::shared_ptr<segmented_AABB> AABB(size_t rows, size_t cols);

    /**
    * Loads a mesh from file. Should by x y z values with no header, space delimited.
    * \param file Fully qualified path to a file.
    */
	void from_json(pt::ptree& mesh);

    /**
    * Writes a mesh to an hdf5 file.
    * \param filename_base Base name of file to write .
    */
	void to_hdf5(std::string filename_base);

    /**
    * Reads a mesh and parameters from an hdf5 file.
    * \param mesh_filename Name of mesh file to read .
    * \param param_filename Name of parameter file to read.
    * \param [OPTIONAL] ic_filename Name of initial condition file to read.
    */
    void from_hdf5(const std::string& mesh_filename,
                   const std::vector<std::string>& param_filename,
                   const std::vector<std::string>& ic_filename);

    /**
     * Reads a partitioned mesh file (.partition) that holds multiple hdf5 files
     * @param partition_filename Name of parition file to read
     * @param cwd Current working directory of the calling code. Because the partition file can either include fully
     *   qualified paths or relative paths, we need to know where we currently are so-as to ensure we open the correct
     *   file.
     */
    void from_partitioned_hdf5(const std::string& partition_filename, boost::filesystem::path cwd = boost::filesystem::path("."));

    /**
    * Sets a new order to the face numbering.
    * \param permutation desired ordering
    */
  void reorder_faces(std::vector<size_t> permutation);

    /**
    * Sets the MPI process ownership of mesh faces and nodes
    */
    void partition_mesh();

    /**
     * Load the partition from the h5 file instead of computing it
     */
    void load_partition_from_mesh(const std::string& mesh_filename);

    /**
    * Figures out which faces lie on the boundary of an MPI process' domain
    */
  void determine_local_boundary_faces();
    /**
    * Figures out which faces are ghosted nearest neighbors of
    * This is useful for the communication needed in nearest-neighbor discretizations of spatial operators.
    */
  void determine_process_ghost_faces_nearest_neighbors();
    /**
    * Set up the comm patterns for communicating nearest neighbor mesh data
    */
  void setup_nearest_neighbor_communication();

  /*
    Debug output for detailed ghost cell_global_id and their owners
  */
  void print_ghost_neighbor_info();

  /**
   * Communicate the variable var for all ghost neighbors
   * This signature supports the use case if _s no-oped to const char * via
   * #define SAFE_CHECKS
   * which is used to debug and print var names. In general the hash variant used via "var_name"_s should
   * be used.
   * @param var Variable name
   */
  void ghost_neighbors_communicate_variable(const std::string& var);
    /**
     * Communicate the variable var for all ghost neighbors
     * @param var Variable name
    */
  void ghost_neighbors_communicate_variable(const uint64_t& var);

  void ghost_to_neighbors_communicate_variable(const std::string& var);
  void ghost_to_neighbors_communicate_variable(const uint64_t& var);

    /**
    * Figures out which faces are required in the ghost region of an MPI process.
    * \param max_distance the maximum distance needed for communication
    */
  void determine_process_ghost_faces_by_distance(double max_distance);
    /**
    * Shrink the local mesh to only contain owned entries and relevant ghost entries
    */
  void shrink_local_mesh_to_owned_and_distance_neighbors();



    /**
     * Serializes a mesh attribute to file so it can be read into the model.
     * \param output_path Full qualified path to a file to output to.
     * \param variable The variable from a module to be output.
     */
    void serialize_parameter(std::string output_path, std::string parameter);
    /**
    * Initializes a triangulation from the given x,y,z vectors
    * \param x X values
    * \param y Y values
    * \param z Z values
    */
//    void init(vector x, vector y, vector z);

    /**
    * Return the number of faces in the local triangluation
    * \return Number of triangle faces
    */
    size_t size_faces();

    /**
    * Return the number of faces in the global triangluation
    * \return Number of triangle faces
    */
    size_t size_global_faces();

    /**
    * Return the number of verticies in the triangulation
    * \return Number of vertices
    */
    size_t size_vertex();

    /**
     * Locates the closest triangle (based on centers) that the point x y lies on in 2 dimensions.
     * This will always return a triangle, even if the query point is outside of the domain.
     * The triangle found does not nessessairly contain the point specified if another triangle's center is closer.
     * @param x x coordinate of query point
     * @param y y coordinate of query point
     * @return
     */
    mesh_elem find_closest_face(double x, double y) const;

    /**
     * Locates the closest triangle (based on centers) that the point x y lies on in 2 dimensions.
     * This will always return a triangle, even if the query point is outside of the domain
     * The triangle found does not nessessairly contain the point specified if another triangle's center is closer.
     * @param query Query point
     * @return
     */
	mesh_elem find_closest_face(Point_2 query) const;

    /**
     * Locates the triangles (based on centers) within a given radius. Example:
     * @code
     *       auto s = global_param->stations().at(0);
     *       auto faces_in_radius = domain->find_faces_in_radius(s->x(),s->y(),1000.);
     *
     *       for (auto& f : faces_in_radius)
     *       {
     *           f->set_face_data("radius_test", 1);
     *       }
     * @endcode
     * @param center x,y-coordinate of the center to search from
     * @param radius radius within which to search. Given in mesh units (e.g., metres)
     * @return Vector of faces
     */
    std::vector< mesh_elem > find_faces_in_radius(Point_2 center, double radius) const;

	/**
	 * Locates the triangles (based on centers) within a given radius.
	 * @param x x-coordinate of the center to search from
	 * @param y y-coordinate of the center to search from
	 * @param radius radius within which to search. Given in mesh units (e.g., metres)
	 * @return Vector of faces
	 */
    std::vector< mesh_elem > find_faces_in_radius(double x, double y, double radius) const;

    /**
     * Lcoates the triangle that contains the query point. Guaranteed that if a triangle is found, the point lies inside the triangle.
     *
     * @param x
     * @param y
     * @return triangle if found, otherwise nullptr
     */
    mesh_elem locate_face(double x, double y);
    mesh_elem locate_face(Point_2 query);

    /**
    * Returns the finite face at index i. A given index will always return the same face.
    * \param i Index
    * \return A face handle to the ith face
    */
    mesh_elem face(size_t i);

    /**
    * Returns the vector of locally owned global IDs
    * - "locally owned" = this rank is responsible for their data
    * \return A vector of global IDs for the locally owned elements
    */
    const std::vector<int>& get_global_IDs() const;

    /**
    * Returns the finite vertex at index i. A given index will always return the same vertex.
    * \param i Index
    * \return A vertex handle to the ith vertex
    */
    Delaunay::Vertex_handle vertex(size_t i);
#ifdef NOMATLAB
    /**
    * If Matlab integration is enabled, plots the given variable at the current timestep.
    * \param ID variable
    */
    void plot(std::string ID);

    /**
    * If Matlab integration is enabled, plots the given variable over time for the given triangle that contains the point x,y
    * \param x x coord
    * \param y y coord
    * \param ID variable to plot
    */

    void plot_time_series(double x, double y, std::string ID);
#endif

    /**
    * Saves the timeseries as a csv for the triangle that contains the point x,y to file
    * \param x x coord
    * \param y y coord
    * \param fname File same
    */
    void timeseries_to_file(double x, double y, std::string fname);

    /**
    * Saves the timeseries as a csv for the given triangle to file
    * \param m mesh element
    * \param fname File same
    */
    void timeseries_to_file(mesh_elem m, std::string fname);

    /**
     * If output to the mesh vtk/vtu format is required, this will be allocate the vtk data structure.
     */
    void init_vtkUnstructured_Grid(std::vector<std::string> output_variables);

    /// Initializes all the face timeseries to hold the selected variables
    /// @param variables
    void init_timeseries(std::set< std::string > variables);

    /// Initializes the face vectors
    /// @param variables
    void init_vectors(std::set<std::string>& variables);

    /// Initializes all the faces to hold module data
    /// @param modules
    void init_module_data(std::set< std::string > modules);

    /// Initalizes all the face-data data structures: variables, module data, vectors.
    /// Can be done individually but this only requires one pass over the triangulation and is thus faster
    /// @param timeseries
    /// @param vectors
    /// @param modules
    /// @param module_data
    void init_face_data(std::set< std::string >& timeseries,
                  std::set< std::string >& vectors,
                  std::set< std::string >& module_data);

    /**
     * Prunes the internal vector that holds faces to only hold a subset. Does not actually remove the faces from the
     * triangulation. Cannot be used with MPI ranks >1 and outside point mode.
     * @param faces
     */
    void prune_faces(std::vector<Face_handle>& faces);

	/**
	 * Updates the internal vtk structure with this timesteps data.
	 * Must be called prior to calling the write_vt* functions.
	 * The write_vt* functions could call this, however it makes them not threadsafe.
	 * If output_variables is empty, it will write all variables out
	 * @param output_variables Selected variables to write out.
	 */
    void update_vtk_data(std::vector<std::string> output_variables);

    /**
    * Saves the mesh with this timesteps values to a vtu file for visualization in Paraview
    */
	void write_vtu(std::string fname);


	/**
	 * Returns true if this is a geogrphic mesh
	 * @return
	 */
	bool is_geographic();

    /**
     * Returns the proj4 description of the projection used
     * @return
     */
	std::string proj4();

	//holds the spatial search tree
	//http://doc.cgal.org/latest/Spatial_searching/index.html
	boost::shared_ptr<Tree> dD_tree;

    /**
     * Set the the private variable for writing parameters in vtu output
     */
    void write_param_to_vtu(bool write_param);
    /**
     * Set the the private variable for writing the ghost neighbor data in vtu output
     */
    void write_ghost_neighbors_to_vtu(bool write_ghost_neighbors);

    /**
     * Returns the set of parameters available on the triangulation
     * @return
     */
    std::set<std::string> parameters();

    bool _terrain_deformed;

    /**
     * Minimium z elevation in mesh
     * @return
     */
    double min_z();

    /**
     * Maximum elevation in mesh
     * @return
     */
    double max_z();

    //Point to the global object that contains paramter information.
    boost::shared_ptr<global> _global;

    //this holds the parameters that we load
    //however, core might have found some parameters from modules
    // it will have to insert them into this list so that the static hashmaps can be properly init
    std::set<std::string> _parameters;

#ifdef USE_MPI
    boost::mpi::environment _mpi_env;
    boost::mpi::communicator _comm_world;
#endif

protected:

    enum GHOST_TYPE
    {
        NONE,
        NEIGH,
        DIST
    };

    class MESH_VERSION
    {
      public:
        MESH_VERSION()
        {
            MAJOR=1;
            MINOR=0;
            PATCH=0;
        }

        int MAJOR;
        int MINOR;
        int PATCH;

        std::string to_string()
        {
            std::string s;
            s = std::to_string(MAJOR) + "." + std::to_string(MINOR) + "." + std::to_string(PATCH);
            return s;
        }
        void from_string(std::string s)
        {
            if( s == "")
            {
                s = "1.0.0";
            }

            std::vector<std::string> r;
            boost::algorithm::split(r, s, boost::is_any_of("."));

            MAJOR = std::stoi(r[0]);
            MINOR = std::stoi(r[1]);
            PATCH = std::stoi(r[2]);
        }


        bool mesh_ver_meets_min_json()
        {
            return MAJOR >=1;
        }

        /**
         * Basic h5 versioning is currently the same as min json
         * @return
         */
        bool mesh_ver_meets_min_h5()
        {
            return MAJOR >=1;
        }
        bool mesh_ver_meets_min_partition()
        {
            return MAJOR >=2 && MINOR >=0;
        }

    } _version;

    /**
     * Loads a given mesh as h5 into the triangulation. Only loads the main topology.
     * @param mesh_filename
     */
    void load_mesh_from_h5(const std::string& mesh_filename);

    void determine_ghost_owners();

    /**
     * Build the spatial search dD tree. Assumes _num_global_faces has been set and this needs to be called after
     * the partition has happened
     */
    void _build_dDtree();

    size_t _num_faces; //number of faces, in MPI mode this will be the local number of faces
    size_t _num_global_faces; //number of global faces
    size_t _num_vertex; //number of rows in the original data matrix. useful for exporting to matlab, etc
    K::Iso_rectangle_2 _bbox;
    bool _is_geographic;
    bool _mesh_is_from_partition;
    int _UTM_zone;

    std::string _srs_wkt;
    //holds the vtk ugrid if we are outputing to vtk formats
    vtkSmartPointer<vtkUnstructuredGrid> _vtk_unstructuredGrid;

	//holds the vectors we use to create the vtu file
#ifdef USE_SPARSEHASH
    google::dense_hash_map< std::string, vtkSmartPointer<vtkFloatArray>  > data;
    google::dense_hash_map< std::string, vtkSmartPointer<vtkFloatArray>  > vectors;
    google::dense_hash_map< std::string, vtkSmartPointer<vtkFloatArray>  > vertex_data;
#else
	std::map<std::string, vtkSmartPointer<vtkFloatArray> > data;
	std::map<std::string, vtkSmartPointer<vtkFloatArray> > vectors;
        std::map<std::string, vtkSmartPointer<vtkFloatArray> > vertex_data;
#endif

    //should we write parameters to the vtu file?
    bool _write_parameters_to_vtu;
    //should we write ghost neighbor faces to the vtu file?
    bool _write_ghost_neighbors_to_vtu;

    // min and max elevations
    double _min_z;
    double _max_z;

    //If the triangulation is traversed using the finite_faces_begin/end iterators, the determinism of the order of traversal is not guaranteed
    //as well, it seems to prevent openmp for applying parallelism to the for-loops. Therefore, we will just store a predefined list of faces and vertex handles
    //that allows us to traverse the triangulation in a deterministic order, as well as play nice with openmp
    std::vector< mesh_elem > _faces;
    std::vector< Delaunay::Vertex_handle > _vertexes;

    // Maps a global index to a local index
    // Size of this vector is the number of locally owned elements in _faces and so incl _local_faces + ghosts
    // key=global_index, entry=local_index
    std::map<int,int> _global_to_locally_owned_index_map;

    // Maps a global index to a local index FOR USE WITH _local_faces
    // Size of this vector is the number of locally owned elements in _local_faces and DOES NOT INCL ghosts
    // key=global_index, entry=local_index
    std::map<int,int> _global_to_local_faces_index_map;

    // All MPI process are aware of the local sizes for all other MPI processes
    // Index = MPI rank
    std::vector<int> _num_faces_in_partition;
    std::vector<int> _local_sizes;

    // If we are not using MPI, some code paths might still want to make use of these
    // This is initialized to [0, num_faces - 1]
    // In MPI mode this contains the global face index start and end
    int  global_cell_start_idx, // in non-MPI this is 0,
        global_cell_end_idx; // in non-MPI this is equal to _num_faces - 1

    // The faces owned by this rank. In non MPI mode, this is identical to _faces
    // does not include any ghosts
    std::vector< mesh_elem > _local_faces;


    std::vector< std::pair<mesh_elem,bool> > _boundary_faces;

    // Array of pointers to the ghost neighbors for this rank
    // These are only the nearest neighbour ghosts. i.e., GHOST_TYPE::NEIGH
    std::vector< mesh_elem > _ghost_neighbors;

    // Array of which other rank owns each ghost neighbor
    std::vector< int > _ghost_neighbor_owners;

    // Array of pointers to all ghost faces (in buffer region) for this rank
    // This is all the ghost faces and is the super-set of ghosts that contains
    // _ghost_neighbors + the distance (type 2) ghosts
    std::vector< mesh_elem > _ghost_faces;

    // The communication partnership for each rank
    // Partner ID, (start_local_idx, length)
    std::map< int, std::pair<int,int> > _comm_partner_ownership;

    // maps a global index to a local ghost
    std::map<int,int> _global_index_to_local_ghost_map; // key=cell_global_id, entry=local index in _ghost_neighbors array

    std::map< int, std::vector<int> > global_indices_to_send; // key=process to send to, entry = vector of global ids to send

  std::map< int, std::vector<int>> local_indices_to_send; // key=process to send to, entry=locally owned index that corresponds to the global index to be sent
  std::map< int, std::vector<mesh_elem> > local_faces_to_send; // key=process to send to, entry=locally owned pointer to face

  std::map< int, std::vector<int>> ghost_indices_to_recv; // key=process to send to, entry=locally owned index that corresponds to the global index to be sent
  std::map< int, std::vector<mesh_elem> >
      ghost_faces_to_recv; // key=process to send to, entry=locally owned pointer to face

    std::vector<int> _global_IDs;

  std::vector< std::shared_ptr<station> > _stations;

  std::string _partition_method;

#ifdef NOMATLAB
    //ptr to the matlab engine
    boost::shared_ptr<maw::matlab_engine> _engine;
    boost::shared_ptr<maw::graphics> _gfx;
#endif

/*
  Datatypes for reading/writing HDF5 files
*/

  // Array datatype for vertices
  hsize_t vertex_dims;
  H5::ArrayType vertex_t;

  // Array datatype for vertices defining faces
  hsize_t elem_dims;
  H5::ArrayType elem_t;

  // Array datatype for neighbors
  hsize_t neighbor_dims;
  H5::ArrayType neighbor_t;

  // Array datatype for proj4 string
  hsize_t proj4_dims;
  H5::StrType proj4_t;

  // Array datatype for mesh version string
  hsize_t meshver_dims;
  H5::StrType meshver_t;

  // Array datatype for partition type string
  hsize_t partition_type_dims;
  H5::StrType partition_type_t;

  // Array datatype for is_geographic
  hsize_t geographic_dims;

  // Array datatype for is_partition
  hsize_t partition_dims;

};

/**
* \typedef mesh
* Provides a convenience typedef for passing around mesh pointers.
*/
typedef boost::shared_ptr<triangulation> mesh;

/**
* \class rect
* Sub-bounding box rectangle.
*/
class rect
{
public:
    /**
    * Creates a bounding box rectangle with the given coordinates
    * \param coord x,y coordinate pairs
    */
	rect( arma::mat* coord )
	{
		this->coord = coord;
	}
	~rect()
	{
		delete coord;
	}
	arma::mat* coord;
    tbb::concurrent_vector<mesh_elem> triangles;
};

/**
* \class segmented_AABB
* \brief segmented axis aligned bounding box
*
* A segmented axis aligned bounding box for quickly spatially binning triangles. A poor man's k-tree.
* 	         4       top      3
*	         +--------+--------+
*	         |        |        |
*	         |        |        |
*	   left  +--------+--------+    right
*	         |        |        |
*	         |        |        |
*	         +-----------------+
*	         1     bottom      2
*	 +y
*	 ^
*	 |
*	 +-> +x
*
*/
class segmented_AABB
{
public:
	segmented_AABB();
	~segmented_AABB();

    /**
    * Creates a segmented AABB for the given triangulation of size rows x cols
    * \param domain A triangulation mesh
    * \param rows number of rows
    * \param cols number of cols
    */
	void make( triangulation* domain,  size_t rows, size_t cols);

    /**
    * Returns the sub rectangle
    * \param row
    * \param col
    * \return The sub rectangle
    */
	rect* get_rect( size_t row, size_t col);

    /**
    * \var n_rows
    * Number of rows in the AABB
    */
	size_t n_rows;

    /**
    * \var n_cols
    * Number of cols in the AABB
    */
	size_t n_cols;

    /**
    * Determins if a point exists within a given rectangle
    * \param v Vertex handle toa  point (vertex) in the triangulation
    * \param r A subrect
    * \returns true if point lies whitin r
    */
    bool pt_in_rect(Delaunay::Vertex_handle v, rect* r);

    /**
    * Determines if a point given by x y coords exists in a given rectangle
    * \param x coord
    * \param y coord
    * \param r sub rectangle
    * \returns true if point is within r
    */
	bool pt_in_rect(double x, double y, rect* r);

private:
	tbb::concurrent_vector<tbb::concurrent_vector<rect*> > m_grid;

};

template < class Gt, class Fb>
std::vector<std::shared_ptr<station>>& face<Gt, Fb>::stations()
{
    return _stations;
}

template < class Gt, class Fb>
std::shared_ptr<station>& face<Gt, Fb>::nearest_station()
{
    return _nearest_station;
}

template < class Gt, class Fb>
bool face<Gt, Fb>::has_parameter(const std::string& variable)
{
    return _parameters.has(variable);
}

template < class Gt, class Fb>
bool face<Gt, Fb>::has_parameter(const uint64_t& hash)
{
    return _parameters.has(hash);

}

template < class Gt, class Fb >
double& face<Gt, Fb>::parameter(const std::string& variable)
{
    return _parameters[variable];
};

template < class Gt, class Fb >
double& face<Gt, Fb>::parameter(const uint64_t& hash)
{
    return _parameters[hash];


};

template < class Gt, class Fb >
bool face<Gt, Fb>::has_initial_condition(std::string key)
{
    return _initial_conditions.has(key);
}

template < class Gt, class Fb >
void face<Gt, Fb>::set_initial_condition(std::string key, double value)
{
    _initial_conditions[key] = value;
}

template < class Gt, class Fb >
double face<Gt, Fb>::get_initial_condition(std::string key)
{
    return _initial_conditions[key];
};

template < class Gt, class Fb >
std::vector<std::string>  face<Gt, Fb>::parameters()
{
    return _parameters.variables();
};

template < class Gt, class Fb >
std::vector<std::string>  face<Gt, Fb>::initial_conditions()
{
    return _initial_conditions.variables();
};

template < class Gt, class Fb >
std::vector<std::string>  face<Gt, Fb>::vectors()
{
    return _module_face_vectors.variables();
};

template < class Gt, class Fb >
face<Gt, Fb>::~face()
{

};
template < class Gt, class Fb >
face<Gt, Fb>::face()
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();
    _center = NULL;
    _normal = NULL;
    _area = -1.;
    _is_geographic = false;



}

template < class Gt, class Fb>
face<Gt, Fb>::face(Vertex_handle v0,
                   Vertex_handle v1,
                   Vertex_handle v2)
        : Fb(v0, v1, v2)
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();
    _center = NULL;
    _normal = NULL;
    _area = -1.;
    _is_geographic = false;

}

template < class Gt, class Fb >
face<Gt, Fb>::face(Vertex_handle v0,
                   Vertex_handle v1,
                   Vertex_handle v2,
                   Face_handle n0,
                   Face_handle n1,
                   Face_handle n2)
        : Fb(v0, v1, v2, n0, n1, n2)
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();
    _center = NULL;
    _normal = NULL;
    _area = -1.;
    _is_geographic = false;

}

template < class Gt, class Fb >
face<Gt, Fb>::face(Vertex_handle v0,
                   Vertex_handle v1,
                   Vertex_handle v2,
                   Face_handle n0,
                   Face_handle n1,
                   Face_handle n2,
                   bool c0,
                   bool c1,
                   bool c2)
        : Fb(v0, v1, v2, n0, n1, n2)
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<timeseries>();
    _center = NULL;
    _normal = NULL;
    _area = -1.;
    _is_geographic = false;


}

template < class Gt, class Fb>
double face<Gt, Fb>::aspect()
{
    if (_azimuth == -1)
    {
        if(!_normal)
            this->normal();

        _azimuth = math::gis::cartesian_to_bearing(Vector_2((*_normal)[0],(*_normal)[1])) * M_PI/180.; //need in radians

    }

    return _azimuth;
}

template < class Gt, class Fb>
Vector_2 face<Gt, Fb>::edge_unit_normal(int i)
{
    auto e = edge(i);
    auto e1 = edge( (i+1) % 3);

    //use this method to get the vector going in the right dir
//    http://gamedev.stackexchange.com/a/26952

    //one of the normals is defined as
    Vector_2 n(e.y(),-e.x());

    //dot normal with e1
    double D = e1.x()*n.x()+e1.y()*n.y();

    //if positive, negate to point the proper dir
    if(D > 0)
        n = -n;

    return n/CGAL::sqrt(n.squared_length());

};

template < class Gt, class Fb>
Vector_2 face<Gt, Fb>::edge(int i)
{
    auto a = this->vertex(_domain->ccw(i))->point();
    auto b = this->vertex(_domain->cw(i))->point();

    auto v = b-a;

    return Vector_2(v.x(),v.y());//,v.z());

};
template < class Gt, class Fb>
std::pair<Point_2,Point_2> face<Gt, Fb>::edge_vertexes(int i)
{
    auto a = this->vertex(_domain->ccw(i))->point();
    auto b = this->vertex(_domain->cw(i))->point();
    return std::make_pair(a,b);
};
template < class Gt, class Fb>
Point_2 face<Gt, Fb>::edge_midpoint(int i)
{
    auto a = this->vertex(_domain->ccw(i))->point();
    auto b = this->vertex(_domain->cw(i))->point();
    Point_2 midpoint( (a.x()+b.x())/2.0 , (a.y() + b.y()) / 2.0);
    return midpoint;
};

template < class Gt, class Fb>
double face<Gt, Fb>::edge_length(int i)
{
    auto e = edge(i);

    return CGAL::sqrt(e.squared_length());
};

template < class Gt, class Fb>
double face<Gt, Fb>::slope()
{
    if (_slope == -1)
    {
        if(!_normal)
            this->normal();

        //z surface normal
        arma::vec n(3);
        n(0) = 0.0; //x
        n(1) = 0.0; //y
        n(2) = 1.0;

        arma::vec normal(3);
        normal(0) = (*_normal)[0];
        normal(1) = (*_normal)[1];
        normal(2) = (*_normal)[2];

        _slope = acos(arma::norm_dot(normal, n));
    }

    return _slope;
}

template < class Gt, class Fb>
bool face<Gt, Fb>::intersects(face<Gt, Fb>::Face_handle fh)
{
    //does t contain my points?
    bool intersect = fh->contains(this->center());
    if (!intersect)
    {
        if (fh->contains(this->vertex(0)->point()) ||
            fh->contains(this->vertex(1)->point()) ||
            fh->contains(this->vertex(2)->point()))
        {
            intersect = true;
        }
    }
    return intersect;
}

template < class Gt, class Fb>
const typename face<Gt, Fb>::Face_handle face<Gt, Fb>::find_closest_face(double azimuth, double distance)
{
    return _domain->find_closest_face(math::gis::point_from_bearing(center(), azimuth, distance));
};

template < class Gt, class Fb>
Vector_3 face<Gt, Fb>::normal()
{
    if(!_normal)
    {
        if(_is_geographic)
        {

//            OGRSpatialReference monUtm;
//
//            OGRSpatialReference monGeo;
//            monGeo.SetWellKnownGeogCS("WGS84");


            CGAL::Point_3<K> v0(this->vertex(0)->point()[0]*100000., this->vertex(0)->point()[1]*100000.,this->vertex(0)->point()[2]);
            CGAL::Point_3<K> v1(this->vertex(1)->point()[0]*100000., this->vertex(1)->point()[1]*100000.,this->vertex(1)->point()[2]);
            CGAL::Point_3<K> v2(this->vertex(2)->point()[0]*100000., this->vertex(2)->point()[1]*100000.,this->vertex(2)->point()[2]);

            _normal = boost::make_shared<Vector_3>(CGAL::unit_normal(v0, v1, v2));

        }
        else
            _normal = boost::make_shared<Vector_3>(CGAL::unit_normal(this->vertex(0)->point(), this->vertex(1)->point(), this->vertex(2)->point()));
    }

    return *_normal;
}

template < class Gt, class Fb>
Point_3 face<Gt, Fb>::center()
{
    if (!_center)
    {
        _center = boost::make_shared<Point_3>(CGAL::centroid(this->vertex(0)->point(), this->vertex(1)->point(), this->vertex(2)->point()));
        _x=_center->x();
        _y=_center->y();
        _z=_center->z();
    }
    //return CGAL::centroid(this->vertex(0)->point(), this->vertex(1)->point(), this->vertex(2)->point());
    return *_center;

}
template < class Gt, class Fb>
bool face<Gt, Fb>::contains(Point_3 p)
{
    return this->contains(p.x(), p.y());
}

template < class Gt, class Fb>
bool face<Gt, Fb>::contains(double x, double y)
{

    double x1 = this->vertex(1)->point().x();
    double y1 = this->vertex(1)->point().y();

    double x2 = this->vertex(2)->point().x();
    double y2 = this->vertex(2)->point().y();

    double x3 = this->vertex(0)->point().x();
    double y3 = this->vertex(0)->point().y();

    double lambda1 = ((y2 - y3)*(x - x3)+(x3 - x2)*(y - y3)) / ((y2 - y3)*(x1 - x3)+(x3 - x2)*(y1 - y3));

    if( !(lambda1 > 0.0 && lambda1 < 1.0) )
        return false; //bail early if possible

    double lambda2 = ((y3 - y1)*(x - x3)+(x1 - x3)*(y - y3)) / ((y3 - y1)*(x2 - x3)+(x1 - x3)*(y2 - y3));
    double lambda3 = 1.0 - lambda1 - lambda2;

    return lambda1 > 0.0 && lambda1 < 1.0
           && lambda2 > 0.0 && lambda2 < 1.0
           && lambda3 > 0.0 && lambda3 < 1.0;

}
template < class Gt, class Fb>
std::vector<std::string> face<Gt, Fb>::variables()
{

    return _variables.variables();
}


template < class Gt, class Fb>
bool face<Gt, Fb>::has(const std::string& variable)
{
    return _variables.has(variable);

};

template < class Gt, class Fb>
bool face<Gt, Fb>::has(const uint64_t& hash)
{
    return _variables.has(hash);
}

template < class Gt, class Fb>
double& face<Gt, Fb>::operator[](const uint64_t& hash)
{
     return _variables[hash];
}

template < class Gt, class Fb>
double& face<Gt, Fb>::operator[](const std::string& variable)
{
    return _variables[variable];
}

template < class Gt, class Fb >
bool  face<Gt, Fb>::has_vegetation()
{
    if (has_parameter("landcover"_s) )
        return true;
    if (has_parameter("canopyType"_s))
        return true;
    if (has_parameter("CanopyHeight"_s))
        return true;

    return false;
}
template < class Gt, class Fb >
double face<Gt, Fb>::veg_attribute(const std::string &variable)
{
    double result = 0;

    // first see if we have a distributed map of this parameter
    if(has_parameter(variable))
    {
        result = parameter(variable);
    }
    else if(has_parameter("landcover"_s)) // Ok, try to look it up in a classified landcover lookup table
    {
        int LC = parameter("landcover"_s);
        auto param = _domain->_global->parameters; //this grabs the loaded landcover map
        result = param.get<double>("landcover." + std::to_string(LC) + "."+variable);
    }
    else
    {
        BOOST_THROW_EXCEPTION(module_error() << errstr_info("Parameter " + variable +" does not exist."));
    }

    return result;
};

template < class Gt, class Fb>
Vector_3 face<Gt, Fb>::face_vector(const std::string& variable)
{
    return _module_face_vectors[variable];
};

template < class Gt, class Fb>
void face<Gt, Fb>::init_time_series(std::set<std::string>& variables)
{
    _variables.init(variables);
}

template < class Gt, class Fb>
void face<Gt, Fb>::init_vectors(std::set<std::string>& variables)
{
    _module_face_vectors.init(variables);
}

template < class Gt, class Fb>
void face<Gt, Fb>::init_parameters(std::set<std::string>& parameters)
{

    _parameters.init(parameters);
}

template < class Gt, class Fb>
void face<Gt, Fb>::init_module_data(std::set<std::string>& modules)
{
    _module_face_data.init(modules);
}

template < class Gt, class Fb>
timeseries::variable_vec face<Gt, Fb>::face_time_series(std::string ID)
{
    return _data->get_time_series(ID);
}

template < class Gt, class Fb>
void face<Gt, Fb>::next()
{
    ++_itr;
}

template < class Gt, class Fb>
void face<Gt, Fb>::reset_to_begining()
{
    _itr = _data->begin();
}

template < class Gt, class Fb>
void face<Gt, Fb>::to_file(std::string fname)
{
    _data->to_file(fname);
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_x()
{
    if(!_center)
        this->center();

    return _x;
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_y()
{
    if(!_center)
        this->center();


    return _y;
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_z()
{
    if(!_center)
        this->center();

    return _z;
}
template < class Gt, class Fb>
boost::shared_ptr<timeseries> face<Gt, Fb>::get_underlying_timeseries()
{
    return _data;
}

template < class Gt, class Fb>
timeseries::iterator face<Gt, Fb>::now()
{
    return _itr;
}


template < class Gt, class Vb>
template<typename T>
T& face<Gt, Vb>::make_module_data(const std::string &module)
{

    //we don't already have this, make a new one.
    if(!_module_face_data[module])
    {
//        T* data = new T;
        _module_face_data[module] = std::make_unique<T>();
    }

    return get_module_data<T&>(module);
}


template < class Gt, class Fb>
template < typename T>
T& face<Gt, Fb>::get_module_data(const std::string &module)
{
    return dynamic_cast<T&>(*_module_face_data[module]);
}

template < class Gt, class Fb>
void face<Gt, Fb>::set_face_vector(const std::string& variable, Vector_3 v)
{
    _module_face_vectors[variable] = v;
};

// I don't think this is used anywhere and is maybe not worth keeping given the make_module_data exists
//template < class Gt, class Fb>
//void face<Gt, Fb>::set_module_data(const std::string &module, face_info *fi)
//{
//    _module_face_data[module] = fi;
//}
template < class Gt, class Fb>
double face<Gt, Fb>::get_area()
{
    if(_area == -1)
    {
        // supports geographic
        if(has_parameter("area"_s))
        {
            _area = parameter("area"_s);
        }
        else
        {

            auto& pa = this->vertex(0)->point();
            auto& pb = this->vertex(1)->point();
            auto& pc = this->vertex(2)->point();


            //same way it's done in mesher for consistency
            typename Fb::Geom_traits traits;
            _area = CGAL::to_double(traits.compute_area_2_object()(pa, pb, pc));

        }
    }


    return _area;
}
template < class Gt, class Fb>
double face<Gt, Fb>::get_subgrid_z(Point_2 query)
{

            auto x1 = this->vertex(0)->point().x();
            auto y1 = this->vertex(0)->point().y();
            auto z1 = this->vertex(0)->point().z();

            auto x2 = this->vertex(1)->point().x();
            auto y2 = this->vertex(1)->point().y();
            auto z2 = this->vertex(1)->point().z();

            auto x3 = this->vertex(2)->point().x();
            auto y3 = this->vertex(2)->point().y();
            auto z3 = this->vertex(2)->point().z();

            double ux = x2-x1;
            double uy = y2-y1;
            double uz = z2-z1;

            double vx = x3-x1;
            double vy = y3-y1;
            double vz = z3-z1;

            double a, b, c;
            a = uy * vz - vy * uz;
            b = vx * uz - ux * vz;
            c = ux * vy - vx * uy;

            //solve for d
            double d =a*x1 + b*y1 + c*z1;

            double z = -(a*query.x()+b*query.y()-d)/c;

            return z;


};

template <typename T>
T determine_owner_of_global_index(T index, std::vector<T> num_faces_in_partition)
{
  assert(num_faces_in_partition.size() != 0 );

  // zero init variables
  T owner = T();
  T highest_owned_index = T();

  for(T i=0; i<static_cast<T>(num_faces_in_partition.size()); ++i)
  {
    highest_owned_index += num_faces_in_partition[i];
    if (index < highest_owned_index)
    {
      owner = i;
      break;
    }
  }
  return owner;
}
