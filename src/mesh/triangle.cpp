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


#include "triangle.h"

size_t triangle::HashCompare::hash(const std::string& x)
{
    boost::crc_32_type crc32;
    std::string xlower = boost::algorithm::to_lower_copy<std::string>(x);
    crc32.process_bytes(xlower.c_str(), xlower.length());

    return crc32.checksum();
}

bool triangle::HashCompare::equal(const std::string& s1, const std::string& s2)
{
    std::string ss1 = boost::algorithm::to_lower_copy<std::string>(s1);
    std::string ss2 = boost::algorithm::to_lower_copy<std::string>(s2);

    return ss1 == ss2;
}

triangle::triangle( point vertex1, point vertex2, point vertex3, size_t cur_rec_depth/*=1*/)
{
	m_cur_rec_depth = cur_rec_depth;
 	m_sub_tri = NULL;
	set_vertex_values(vertex1,vertex2, vertex3);

	radiation_dir = 0.0;
	radiation_diff = 0.0;
	shadow = 0.0;
	z_prime = 0.0;

	m_slope = 0.0;
	m_azimuth = 0.0;

}

triangle::triangle(size_t cur_rec_depth)
{
	m_sub_tri = NULL;

	m_cur_rec_depth = cur_rec_depth;

	radiation_dir = 0.0;
	radiation_diff = 0.0;
	shadow = 0.0;
	z_prime = 0.0;

	m_slope = 0.0;
	m_azimuth = 0.0;
}


bool triangle::contains( point xy )
{
	return contains(xy.x,xy.y);
}

double triangle::get_x()
{
    return center.x;
}
double triangle::get_y()
{
    return center.y;
}
double triangle::get_z()
{
    return center.z;
}     
boost::posix_time::ptime triangle::get_ptime()
{
    return _current_time;
}
void triangle::set_current_time(boost::posix_time::ptime time)
{
    _current_time = time;
}
void triangle::set_vertex_values( point vertex1, point vertex2, point vertex3)
{
	m_vertex_list[0].x = vertex1.x;
	m_vertex_list[0].y = vertex1.y;
	m_vertex_list[0].z = vertex1.z;

	m_vertex_list[1].x = vertex2.x;
	m_vertex_list[1].y = vertex2.y;
	m_vertex_list[1].z = vertex2.z;

	m_vertex_list[2].x = vertex3.x;
	m_vertex_list[2].y = vertex3.y;
	m_vertex_list[2].z = vertex3.z;

	//center of the triangle
	point pos = calc_center(this);		

	m_center.x = pos.x;
	m_center.y = pos.y;
// 	std::cout.precision(16);
// 	std::cout << "vertex" <<std::endl;
// 	std::cout << m_vertex_list[0].x <<" " << m_vertex_list[0].y << std::endl
// 			  << m_vertex_list[1].x <<" " << m_vertex_list[1].y << std::endl
// 			  << m_vertex_list[2].x <<" " << m_vertex_list[2].y << std::endl;
// 	std::cout << "center" <<std::endl;
// 	std::cout << pos.x << " " << pos.y << std::endl;

	if (m_cur_rec_depth != 0)
	{
		update_subtri();
	}
	
}

void triangle::add_face_data(std::string ID, double data)
{
    face_data::accessor a;

    if(!_data.insert(a,ID)) //insert creates a new pair if not found
    {
        BOOST_THROW_EXCEPTION(mesh_insertion_error()
                    << errstr_info(std::string("Failed to add face data: ") + ID)
                    );
    }

    a->second.push_back(data);
    
//    LOG_DEBUG << a->second[0];
}

double triangle::get_face_data(std::string ID)
{
    face_data::accessor a;
    if(!_data.find(a,ID)) 
    {
        BOOST_THROW_EXCEPTION(mesh_lookup_error()
                    << errstr_info(std::string("No face data: ") + ID)
                    );
    }
    
    double data = a->second[0];//*(a->second.end());
    return  data;

}

void triangle::update_subtri()
{
	int longest;
 	if (m_sub_tri)
	{
		for(int i = 0;i<4;i++)
		{
			delete m_sub_tri[i];
		}
 		delete[] m_sub_tri;	
	}

	// set each sub
	m_sub_tri = new triangle*[4];
	for(int i = 0; i<4;i++)
	{
		m_sub_tri[i] = new triangle(m_cur_rec_depth-1);
	}

// 	std::cout.precision(16);
// 	std::cout << "Main triangle" << std::endl;
// 	std::cout << "[" << m_vertex_list[0].x << " " << m_vertex_list[0].y << "];" << std::endl
// 				<< "[" << m_vertex_list[1].x << " " << m_vertex_list[1].y << "];" << std::endl
// 				<< "[" << m_vertex_list[2].x << " " << m_vertex_list[2].y << "];" << std::endl;

	arma::vec l_sides(3);
	//0-1
	l_sides[0] = sqrt( ((m_vertex_list[0].x)-(m_vertex_list[1].x))*((m_vertex_list[0].x)-(m_vertex_list[1].x)) + (((m_vertex_list[0].y)-(m_vertex_list[1].y))*((m_vertex_list[0].y)-(m_vertex_list[1].y))) * 1.0);
	//1-2
	l_sides[1] = sqrt( ((m_vertex_list[1].x)-(m_vertex_list[2].x))*((m_vertex_list[1].x)-(m_vertex_list[2].x)) + (((m_vertex_list[1].y)-(m_vertex_list[2].y))*((m_vertex_list[1].y)-(m_vertex_list[2].y)))* 1.0);
	//0-2
	l_sides[2] = sqrt( ((m_vertex_list[0].x)-(m_vertex_list[2].x))*((m_vertex_list[0].x)-(m_vertex_list[2].x)) + (((m_vertex_list[0].y)-(m_vertex_list[2].y))*((m_vertex_list[0].y)-(m_vertex_list[2].y)))* 1.0);

	if(l_sides[0] > l_sides[1] || l_sides[0] > l_sides[2])
	{
		longest = 0;
		//midpoint of the longest edge (which is 0-1)
		point *midpt_01 = midpoint(m_vertex_list[0],m_vertex_list[1]);
		
		//midpoint of the opposite edge (which is 0-2)
		point* midpt_02 =midpoint(m_vertex_list[0],m_vertex_list[2]);


		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = midpoint(m_vertex_list[1],m_vertex_list[2]);
	
		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02-01
		//t2 = 01-2-02-01
		//t3 = 01-12-2-01
		//t4 = 01-1-12-01

		//t1 = 0-01-02-01
		point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		point v2;
		v2.x = (midpt_01->x);
		v2.y = (midpt_01->y);	
		point v3;
		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01-2-02-01
		//ptr_point v1;
		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);
		//ptr_point v2;
		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;
		//ptr_point v3;
		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 01-12-2-01
		//ptr_point v1;
		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);
		//ptr_point v2;
		v2.x = (midpt_12->x);
		v2.y = (midpt_12->y);
		//ptr_point v3;
		v3.x = m_vertex_list[2].x;
		v3.y = m_vertex_list[2].y;
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 01-1-12-01
		//v1;
		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);
		//ptr_point v2;
		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;
		//ptr_point v3;
		v3.x = (midpt_12->x);
		v3.y = (midpt_12->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);

		

	}
	else if(l_sides[1] > l_sides[2] || l_sides[1] > l_sides[0])
	{
		longest = 1;

		//midpoint of the longest edge (which is 1-2)
		point* midpt_12 = midpoint(m_vertex_list[1],m_vertex_list[2]);


		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = midpoint(m_vertex_list[0],m_vertex_list[1]);

		//midpoint of the opposite edge2 (which is 0-2)
		point* midpt_02 = midpoint(m_vertex_list[0],m_vertex_list[2]);

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-12-0
		//t2 = 01-1-12-01
		//t3 = 12-2-02-12
		//t4 = 0-12-02-0

		//t1 = 0-01-12-0
		point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		point v2;
		v2.x = (midpt_01->x);
		v2.y = (midpt_01->y);	
		point v3;
		v3.x = (midpt_12->x);
		v3.y = (midpt_12->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01-1-12-01

		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);

		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;	

		v3.x = (midpt_12->x);
		v3.y = (midpt_12->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 12-2-02-12

		v1.x = (midpt_12->x);
		v1.y = (midpt_12->y);

		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;	

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 0-12-02-0

		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;

		v2.x = (midpt_12->x);
		v2.y = (midpt_12->y);	

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);


	}
	else 
	{
		longest = 2;

		//midpoint of the longest edge (which is 0-2)
		point* midpt_02 = midpoint(m_vertex_list[0],m_vertex_list[2]);

		//midpoint of the opposite edge (which is 0-1)
		point* midpt_01 = midpoint(m_vertex_list[0],m_vertex_list[1]);

		//midpoint of the opposite edge2 (which is 1-2)
		point* midpt_12 = midpoint(m_vertex_list[1],m_vertex_list[2]);

		//have 4 sub triangles now:
		//define CCW, left->right
		//t1 = 0-01-02
		//t2 = 01 -1 -02
		//t3 = 1 - 12-02
		//t4 = 12 - 2 -02

		//t1 =0-01-02
		point v1;
		v1.x = m_vertex_list[0].x;
		v1.y = m_vertex_list[0].y;
		point v2;
		v2.x = (midpt_01->x);
		v2.y = (midpt_01->y);	
		point v3;
		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[0]->set_vertex_values(v1,v2,v3);

		//t2 = 01 -1 -02

		v1.x = (midpt_01->x);
		v1.y = (midpt_01->y);

		v2.x = m_vertex_list[1].x;
		v2.y = m_vertex_list[1].y;

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[1]->set_vertex_values(v1,v2,v3);

		//t3 = 1 - 12-02

		v1.x = m_vertex_list[1].x;
		v1.y = m_vertex_list[1].y;

		v2.x = (midpt_12->x);
		v2.y = (midpt_12->y);

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[2]->set_vertex_values(v1,v2,v3);

		//t4 = 12 - 2 -02
		v1.x = (midpt_12->x);
		v1.y = (midpt_12->y);

		v2.x = m_vertex_list[2].x;
		v2.y = m_vertex_list[2].y;

		v3.x = (midpt_02->x);
		v3.y = (midpt_02->y);
		m_sub_tri[3]->set_vertex_values(v1,v2,v3);

	}	
}
triangle& triangle::sub_tri(size_t t)
{
	return *(m_sub_tri[t]);
}


int triangle::intersects( triangle* t )
{
	//I have no children
	if(m_sub_tri == NULL)
	{
		//does t contain my points?
		bool intersect = t->contains(this->get_center());
		if(!intersect)
		{
			if( t->contains(this->get_vertex(0)) ||
				t->contains(this->get_vertex(1)) ||
				t->contains(this->get_vertex(2)))
			{
				intersect = true;
			}
		}
		return intersect == true ? 1:0;
		
	}
	else
		//i have children
	{
		int sum=0;
		for(size_t i = 0; i<4; i++)
		{
			sum += m_sub_tri[i]->intersects(t);
		}
		return sum;
	}
}

bool triangle::contains(double x, double y)
{

	double x1=m_vertex_list[0].x;
	double y1=m_vertex_list[0].y;

	double x2=m_vertex_list[1].x;
	double y2=m_vertex_list[1].y;

	double x3=m_vertex_list[2].x;
	double y3=m_vertex_list[2].y;


	double lambda1= ((y2-y3)*(x-x3)+(x3-x2)*(y-y3))/((y2-y3)*(x1-x3)+(x3-x2)*(y1-y3));
	double lambda2= ((y3-y1)*(x-x3)+(x1-x3)*(y-y3))/((y3-y1)*(x2-x3)+(x1-x3)*(y2-y3));
	double lambda3= 1.0-lambda1-lambda2;

	return lambda1 > 0.0 && lambda1 < 1.0 
		&& lambda2 > 0.0 && lambda2 < 1.0
		&& lambda3 > 0.0 && lambda3 < 1.0;

}

point triangle::get_vertex( size_t vertex )
{
	return m_vertex_list[vertex];
}
// 
// point triangle::operator()(size_t v )
// {
// 	return m_vertex_list[vertex];
// }


point triangle::get_center()
{
	return m_center;
}

point* triangle::midpoint( point& p1, point& p2 )
{
	point* mp = new point;
	//new x value 1/2 way there
	mp->x = p1.x + (p2.x - p1.x)/2.0;
	//slope
	double m = (p2.y-p1.y)/(p2.x-p1.x);
	double b = p1.y-m*p1.x;
	mp->y = m*mp->x+b;

	return mp;
}

point triangle::calc_center( triangle* t )
{

	arma::vec Pa(3);
	Pa(0) = (t->m_vertex_list[0].x);
	Pa(1) = (t->m_vertex_list[0].y);
	Pa(2) = 0;

	arma::vec Pb(3);
	Pb(0) = (t->m_vertex_list[1].x);
	Pb(1) = (t->m_vertex_list[1].y);
	Pb(2) = 0;

	arma::vec Pc(3);
	Pc(0) = (t->m_vertex_list[2].x);
	Pc(1) = (t->m_vertex_list[2].y);
	Pc(2) = 0;

	arma::vec AB = Pb - Pa;
	arma::vec AC = Pc - Pa;
	arma::vec BC = Pc - Pb;

//circumcenter
// 	arma::vec N = arma::cross(AC,AB);
// 	arma::vec L1 = arma::cross(AB,N);
// 	arma::vec L2 = arma::cross(BC,N);
// 	arma::vec P21 = (Pc - Pa)/2;
// 	arma::vec P1 = (Pa + Pb)/2;

	
//incenter

	arma::vec uab  = AB / arma::norm(AB,1);
	arma::vec uac  = AC / arma::norm(AC,1);
	arma::vec ubc  = BC / arma::norm(BC,1);
	arma::vec uba = -uab;

	arma::vec L1 = uab + uac;
	arma::vec L2 = uba + ubc;
	arma::vec P21 = Pb-Pa;
	arma::vec P1 = Pa;
     
	arma::mat ML1(L1);
	arma::mat ML2(L2);

	arma::mat ML = arma::join_rows(ML1,-ML2);

	arma::vec lambda = arma::solve(ML,P21);

	arma::vec pos = P1+lambda(0)*L1;

	point p;
	p.x = pos(0);
	p.y = pos(1);
	return p;
}

void triangle::set_facenormal( arma::vec& normal )
{
	m_surface_normal = normal;
}

arma::vec triangle::get_facenormal()
{
	return m_surface_normal;
}

void triangle::compute_azimuth()
{
	//convert normal to spherical
	double r = sqrt(m_surface_normal(0)*m_surface_normal(0) + m_surface_normal(1)*m_surface_normal(1) + m_surface_normal(2)*m_surface_normal(2));
	double theta = acos(m_surface_normal(2)/r); 

	//y=north
	double phi = atan2(m_surface_normal(1),m_surface_normal(0)); // + M_PI /*-3.14159/2*/; //south == 0
	m_azimuth = phi - M_PI/2.0; //set north = 0

	if(m_azimuth < 0.0)
		m_azimuth += 2.0*M_PI;
}

void triangle::compute_slope()
{
	//z surface normal
	arma::vec n(3);
	n(0) = 0.0; //x
	n(1) = 0.0; //y
	n(2) = 1.0;

	m_slope = acos( arma::norm_dot(m_surface_normal,n));
}

double triangle::azimuth()
{
	return m_azimuth;
}

double triangle::slope()
{
	return m_slope;
}