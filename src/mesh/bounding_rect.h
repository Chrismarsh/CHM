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

#include "triangle.h"
#include <libmaw.h>

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
	std::vector<triangle*> triangles;
	std::vector<int> m_globalID;
};
class bounding_rect
{
public:
	bounding_rect(maw::matlab_engine* m_engine);
	~bounding_rect();
	void make(const maw::d_vec x, const maw::d_vec y,  int n_rows, int n_cols);
	rect* bounding_rect::get_rect( int i, int j );
	int n_rows;
	int n_cols;
	bool pt_in_rect(double x, double y, rect* r);

private:
	//std::vector<rect*> m_rectangles;
	std::vector<std::vector<rect*> > m_grid;
	maw::matlab_engine* m_engine;
};


