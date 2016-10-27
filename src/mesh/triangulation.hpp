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

//for valgrind, remove
#define CGAL_DISABLE_ROUNDING_MATH_CHECK

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <set>
#include <stack>
#include <fstream>
#include <utility>
#define ARMA_DONT_USE_CXX11 //intel on linux breaks otherwise
#include <armadillo>

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>

//#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Triangulation_2.h>

//#define CGAL_LAPACK_ENABLED
//#define CGAL_EIGEN3_ENABLED
//#include <CGAL/Monge_via_jet_fitting.h>
#include <tbb/concurrent_vector.h>

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

//http://www.paraview.org/Bug/print_bug_page.php?bug_id=14164
//http://review.source.kitware.com/#/c/11956/
//until 6.0.1 comes out
//#define VTK_HAS_STD_ISNAN
//#define VTK_HAS_STD_ISINF

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h> //pointdata
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkProbeFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkPlaneSource.h>
#include <vtkGeometryFilter.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkElevationFilter.h>
#include <vtkCurvatures.h>
#include <vtkXMLPolyDataWriter.h>
#ifdef NOMATLAB
#include "libmaw.h"
#endif

#include "vertex.hpp"
#include "face.hpp"

#include "math/distance.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Vector_3 Vector_3;
typedef CGAL::Projection_traits_xy_3<K> Gt; //allows for using 2D algorithms on the 3D points

typedef ex_vertex<Gt> Vb; //custom vertex class
typedef face<Gt> Fb; //custom face class
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Delaunay;

typedef Delaunay::Face_handle mesh_elem;
typedef boost::shared_ptr<tbb::concurrent_vector<double>  > vector;

//search tree typedefs
typedef K::Point_2 Point_2;
typedef boost::tuple<Point_2, Delaunay::Face_handle > Point_and_face;
typedef CGAL::Search_traits_2<K>                       Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_face,
		CGAL::Nth_of_tuple_property_map<0, Point_and_face>,
		Traits_base>                                              Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>          K_neighbor_search;
typedef K_neighbor_search::Tree                             Tree;
typedef K_neighbor_search::Distance                         Distance;

//fwd decl
class segmented_AABB;

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
	std::set<std::string> from_json(pt::ptree& mesh);

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
    * Return the number of faces in the triangluation
    * \return Number of triangle faces
    */
    size_t size_faces();

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
    mesh_elem find_closest_face(double x, double y);

    /**
     * Locates the closest triangle (based on centers) that the point x y lies on in 2 dimensions.
     * This will always return a triangle, even if the query point is outside of the domain
     * The triangle found does not nessessairly contain the point specified if another triangle's center is closer.
     * @param query Query point
     * @return
     */
	mesh_elem find_closest_face(Point_2 query);

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
    Delaunay::Face_handle face(size_t i);

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

    void write_vtp(std::string file_name);

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
private:
    size_t _num_faces; //number of faces
    size_t _num_vertex; //number of rows in the original data matrix. useful for exporting to matlab, etc
    K::Iso_rectangle_2 _bbox;
	bool _is_geographic;
	int _UTM_zone;
	std::string _srs_wkt;
	//holds the vtk ugrid if we are outputing to vtk formats
	vtkSmartPointer<vtkUnstructuredGrid> _vtk_unstructuredGrid;

	//holds the vectors we use to create the vtu file
	std::map<std::string, vtkSmartPointer<vtkFloatArray> > data;
	std::map<std::string, vtkSmartPointer<vtkFloatArray> > vectors;


    //If the triangulation is traversed using the finite_faces_begin/end iterators, the determinism of the order of traversal is not guaranteed
    //as well, it seems to prevent openmp for applying parallism to the for-loops. Therefore, we will just store a predefined list of faces and vertex handles
    //that allows us to traverse the triangulation in a deterministic order, as well as play nice with openmp
    std::vector< Delaunay::Face_handle  > _faces;
    std::vector< Delaunay::Vertex_handle > _vertexes;
    
#ifdef NOMATLAB
    //ptr to the matlab engine
    boost::shared_ptr<maw::matlab_engine> _engine;
    boost::shared_ptr<maw::graphics> _gfx;
#endif

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

