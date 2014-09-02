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


#ifdef NOMATLAB
#include "libmaw.h"
#endif

#include "vertex.hpp"
#include "face.hpp"
//#include "AABB.h"


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;


typedef CGAL::Projection_traits_xy_3<K> Gt; //allows for using 2D algorithms on the 3D points

typedef ex_vertex<Gt> Vb; //custom vertex class
typedef face<Gt> Fb; //custom face class
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds; //our data structure that is using the custom classes
typedef CGAL::Delaunay_triangulation_2<Gt, Tds> Delaunay; //specify a delauany triangulation

typedef Delaunay::Face_handle mesh_elem;
typedef boost::shared_ptr<tbb::concurrent_vector<double> > vector;

//fwd decl
class segmented_AABB;

class triangulation
: public Delaunay
{
public:
    triangulation();
#ifdef MATLAB
    triangulation(boost::shared_ptr<maw::matlab_engine> engine);
#endif
    ~triangulation();

    boost::shared_ptr<segmented_AABB> AABB(size_t rows, size_t cols);
    
    void from_file(std::string file);
    void init(vector x, vector y, vector z);

    //return the number of faces in the triangluation
    size_t size_faces();
    size_t size_vertex();
    mesh_elem locate_face(double x, double y);
    Delaunay::Face_handle face(size_t i);
    Delaunay::Vertex_handle vertex(size_t i);
#ifdef NOMATLAB
    void plot(std::string ID);
    void plot_time_series(double x, double y, std::string ID);
#endif

    void to_file(double x, double y, std::string fname);
    void to_file(mesh_elem m, std::string fname);
    void to_vtu(std::string fname);
private:
    size_t _num_faces; //number of faces
    size_t _num_vertex; //number of rows in the original data matrix. useful for exporting to matlab, etc
    K::Iso_rectangle_2 _bbox;
    
    tbb::concurrent_vector< Delaunay::Face_handle > _faces;
    tbb::concurrent_vector< Delaunay::Vertex_handle > _vertexes;
    
#ifdef NOMATLAB
    //ptr to the matlab engine
    boost::shared_ptr<maw::matlab_engine> _engine;
    boost::shared_ptr<maw::graphics> _gfx;
#endif

};

typedef boost::shared_ptr<triangulation> mesh;


class rect
{
public:
	rect( arma::mat* coord )
	{
		this->coord = coord;
	}
	~rect()
	{
		delete coord;
	}
	arma::mat* coord;
	std::vector<mesh_elem> triangles;
};

// 
// 
// 	     -----------------------
// 
// 	         4       top      3
// 	         +--------+--------+
// 	         |        |        |
// 	         |        |        |
// 	   left  +--------+--------+    right
// 	         |        |        |
// 	         |        |        |
// 	         +-----------------+
// 	         1     bottom      2
// 	 +y
// 	 ^
// 	 |
// 	 +-> +x


class segmented_AABB
{
public:
	segmented_AABB();
	~segmented_AABB();
	void make( triangulation* domain,  size_t rows, size_t cols);
	rect* get_rect( size_t i, size_t j );
	int n_rows;
	int n_cols;
        bool pt_in_rect(Delaunay::Vertex_handle v, rect* r);
	bool pt_in_rect(double x, double y, rect* r);

private:
	std::vector<std::vector<rect*> > m_grid;

};

