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

#define ARMA_DONT_USE_CXX11 //intel on linux breaks otherwise
#include <armadillo>

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/bounding_box.h>
#include <tbb/concurrent_vector.h>

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
#ifdef NOMATLAB
#include "libmaw.h"
#endif

#include "vertex.hpp"
#include "face.hpp"



typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;


typedef CGAL::Projection_traits_xy_3<K> Gt; //allows for using 2D algorithms on the 3D points

typedef ex_vertex<Gt> Vb; //custom vertex class
typedef face<Gt> Fb; //custom face class
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds; //our data structure that is using the custom classes
typedef CGAL::Delaunay_triangulation_2<Gt, Tds> Delaunay; //specify a delauany triangulation

typedef Delaunay::Face_handle mesh_elem;
typedef boost::shared_ptr<tbb::concurrent_vector<double>  > vector;

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
    void from_file(std::string file);

    /**
    * Initializes a triangulation from the given x,y,z vectors
    * \param x X values
    * \param y Y values
    * \param z Z values
    */
    void init(vector x, vector y, vector z);

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
    * Locates the triangle that the point x y lies on in 2 dimensions.
    * \returns NULL if a containing triangle is not found, otherwise returns the triangle
    */
    mesh_elem locate_face(double x, double y);

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

    vtkSmartPointer<vtkUnstructuredGrid> mesh_to_vtkUstructuredGrid();

    /**
    * Saves the mesh with this timesteps values to a vtu file for visualization in Paraview
    */
    void mesh_to_vtu(std::string fname);

    void mesh_to_ascii(std::string file_name);
private:
    size_t _num_faces; //number of faces
    size_t _num_vertex; //number of rows in the original data matrix. useful for exporting to matlab, etc
    K::Iso_rectangle_2 _bbox;

    //If the triangulation is traversed using the finite_faces_begin/end iterators, the determinism of the order of traversal is not guaranteed
    //as well, it seems to prevent openmp for applying parallism to the for-loops. Therefore, we will just store a predefined list of faces and vertex handles
    //that allows us to traverse the triangulation in a deterministic order, as well as play nice with openmp
    tbb::concurrent_vector< Delaunay::Face_handle  > _faces;
    tbb::concurrent_vector< Delaunay::Vertex_handle > _vertexes;
    
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

