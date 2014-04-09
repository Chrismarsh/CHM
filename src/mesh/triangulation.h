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

#include <vector>
#include <iostream>
#include <armadillo>

#include "libmaw.h"
#include "triangle.h"

class triangulation
{

public:
	triangulation(maw::matlab_engine* engine);
	~triangulation();

	//Create a delaunay triangulation
	void create_delaunay(arma::vec* x, arma::vec* y, arma::vec* z);

	//return the size of the triangluation
	size_t size();

	//set the vertex data. It is assumed that t-th triangle's global-id is an index into data
	void set_vertex_data(arma::mat& data);

	//returns the t-th triangle
	triangle& operator()(size_t t);
	
	//computer the face normals for all the triangles.
	void compute_face_normals();

	triangle* find_containing_triangle(double x, double y);

	//creates the triangulation matrix for matlab
	maw::d_mat mf_tri_matrix();
        maw::d_mat mf_elevation_data();
        maw::d_vec mf_face_data(std::string ID);
        

private:
	//the triangles
	std::vector<triangle*> m_triangles;

	size_t m_size;  //number of triangulations
        size_t m_data_size;

	//ptr to the matlab engine
	maw::matlab_engine* m_engine;
        
        std::string ID;

};




