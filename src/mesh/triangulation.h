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


#include <armadillo>

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/bounding_box.h>
#include <tbb/concurrent_vector.h>

#include "libmaw.h"
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
typedef boost::shared_ptr<tbb::concurrent_vector<double> > vector;


class triangulation
: public Delaunay
{
public:
    triangulation();
    triangulation(boost::shared_ptr<maw::matlab_engine> engine);
    ~triangulation();

    void from_file(std::string file);
    void init(vector x, vector y, vector z);

    //return the size of the triangluation
    size_t size();

    mesh_elem locate_face(double x, double y);

    void plot(std::string ID);
    void plot_time_series(double x, double y, std::string ID);
    
    void to_file(double x, double y, std::string fname);
    void to_file(mesh_elem m, std::string fname);
private:
    size_t _size; //number of faces
    size_t _data_size; //number of rows in the original data matrix. useful for exporting to matlab,e tc
    K::Iso_rectangle_2 _bbox;
    //ptr to the matlab engine
    boost::shared_ptr<maw::matlab_engine> _engine;
    boost::shared_ptr<maw::graphics> _gfx;

};

typedef boost::shared_ptr<triangulation> mesh;



