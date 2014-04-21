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


#include <armadillo>


#include <vector>
#include <iostream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>

#include <tbb/concurrent_hash_map.h>
#include <tbb/concurrent_vector.h>


#include <boost/crc.hpp>      // for boost::crc_basic, boost::crc_optimal
#include <boost/cstdint.hpp>  // for boost::uint16_t
#include <boost/utility.hpp>

#include <boost/date_time/posix_time/posix_time.hpp>


#include "point.h"
#include "exception.hpp"
#include "logger.h"
#include "crc_hash_compare.hpp"



class triangle
{
    
private:

        
        typedef tbb::concurrent_vector< double > ts;
        typedef tbb::concurrent_hash_map<std::string, ts ,crc_hash_compare> face_data;
        
	//list of the vertexes
	point m_vertex_list[3];

	//center of the triangle
	point m_center;

	//always 4 subtriangles at the moment
	triangle** m_sub_tri;

	//set the number of sub triangles
	size_t m_cur_rec_depth;

	//surface normal vector
	arma::vec m_surface_normal;

	//helper functions:

	//get the mid point of a line segment
	point* midpoint(point& p1, point& p2);
	point  calc_center(triangle* t);


	//in radians.
	double m_slope;

	//positive, clockwise, from north
	double m_azimuth;
        
        boost::posix_time::ptime _current_time;


        face_data _data;
        
public:
	//use xyz triples in the vector
	//store the index like matlab
	triangle(point vertex1, point vertex2, point vertex3, size_t cur_rec_depth=1);
	triangle(size_t cur_rec_depth);


	/*Setters*/

	//recomputes the subtriangles
	void update_subtri();
	void set_vertex_values( point vertex1, point vertex2, point vertex3);
	void set_facenormal(arma::vec& normal);
        void set_current_time(boost::posix_time::ptime time);
	/*Getters*/

	point get_vertex(size_t vertex);
	point get_center();
        double get_x();
        double get_y();
        double get_z();

        
	double azimuth();
	double slope();

	void compute_azimuth();
	void compute_slope();

	//get the t-th subtriangle
	triangle& sub_tri(size_t t);


	bool contains(double x, double y);
	bool contains(point xy);
	int intersects(triangle* t);
	arma::vec get_facenormal();

        void add_face_data(std::string ID, double data);
        double get_face_data(std::string ID);
        boost::posix_time::ptime get_ptime();
	//information for the physical model
	double radiation_dir;
	double radiation_diff;
	double shadow;
	double z_prime;
	double area;
	double cosi;

	point center;
	point rot_center;
	//this is the index used by matlab's triangulation
	//it is [1 .. N] where N is number of triangles
	//it starts at 1 because Matlab's indexing starts at 1
	size_t global_id[3];

	//we are some i-th triangle.
	int    triangle_id;
		


};

typedef triangle mesh_elem;