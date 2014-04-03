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

class point
{
public:
	point()
	{
		x=0;
		y=0;
		z=0;
	}
	point(double x,double y,double z)
	{
		this->x=x;
		this->y=y;
		this->z=z;
	}

	double x;
	double y;
	double z;
};

class ptr_point
{
public:
	ptr_point()
	{
		x=0;
		y=0;
		z=0;
	}

	double* x;
	double* y;
	double* z;
};

