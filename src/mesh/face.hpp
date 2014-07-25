#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_ds_face_base_2.h>

#include <CGAL/Triangulation_face_base_2.h>

#include <set>
#include <string>

#include <boost/shared_ptr.hpp>

#include "timeseries.hpp"
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_3 Point_3;
    typedef K::Vector_3 Vector_3;

template < class Gt, class Fb = CGAL::Triangulation_face_base_2<Gt> >
class face
: public Fb
{
public:

    typedef typename Fb::Vertex_handle Vertex_handle;
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

    double azimuth();
    double slope();
    Vector_3 normal();
    Point_3 center();
    bool contains(double x, double y);
    void set_face_data(std::string var_ID, double data);
    double face_data(std::string var_ID);
    void init_time_series(std::set<std::string> variables, time_series::date_vec datetime, int size);
    time_series::variable_vec face_time_series(std::string ID);
    void next();
    void reset_to_begining();
    double get_x();
    double get_y();
    double get_z();

private:

    double _slope;
    double _azimuth;
    Vector_3 _normal;
    Point_3 _center;

    boost::shared_ptr<time_series> _data;
    time_series::iterator _itr;


};

template < class Gt, class Fb >
face<Gt, Fb>::face()
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<time_series>();
}

template < class Gt, class Fb>
face<Gt, Fb>::face(Vertex_handle v0,
        Vertex_handle v1,
        Vertex_handle v2)
: Fb(v0, v1, v2)
{
    _slope = -1;
    _azimuth = -1;
    _data = boost::make_shared<time_series>();
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
    _data = boost::make_shared<time_series>();
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
    _data = boost::make_shared<time_series>();
}

template < class Gt, class Fb>
double face<Gt, Fb>::azimuth()
{
    if (_azimuth == -1)
    {
        //convert normal to spherical
	double r = sqrt(_normal[0]*_normal[0] + _normal[1]*_normal[1] + _normal[2]*_normal[2]);
	double theta = acos(_normal[2]/r); 

	//y=north
	double phi = atan2(_normal[1],_normal[0]); // + M_PI /*-3.14159/2*/; //south == 0
	_azimuth = phi - M_PI/2.0; //set north = 0

	if(_azimuth < 0.0)
		_azimuth += 2.0*M_PI;
    }

    return _azimuth;
}

template < class Gt, class Fb>
double face<Gt, Fb>::slope()
{
    if (_slope == -1)
    {
        //z surface normal
	arma::vec n(3);
	n(0) = 0.0; //x
	n(1) = 0.0; //y
	n(2) = 1.0;
        
        arma::vec normal(3);
        normal(0) = _normal[0];
        normal(1) = _normal[1];
        normal(2) = _normal[2];

	_slope = acos( arma::norm_dot(normal,n));
    }

    return _slope;
}

template < class Gt, class Fb>
Vector_3 face<Gt, Fb>::normal()
{
    return _normal;
}

template < class Gt, class Fb>
Point_3 face<Gt, Fb>::center()
{
    return _center;
}

template < class Gt, class Fb>
bool face<Gt, Fb>::contains(double x, double y)
{
    return false;
}

template < class Gt, class Fb>
void face<Gt, Fb>::set_face_data(std::string var_ID, double data)
{
    _itr->set(var_ID, data);
}

template < class Gt, class Fb>
double face<Gt, Fb>::face_data(std::string var_ID)
{
    return _itr->get(var_ID);
}

template < class Gt, class Fb>
void face<Gt, Fb>::init_time_series(std::set<std::string> variables, time_series::date_vec datetime, int size)
{
    _data->init(variables, datetime, size);
    _itr = _data->begin();
}

template < class Gt, class Fb>
time_series::variable_vec face<Gt, Fb>::face_time_series(std::string ID)
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
double face<Gt, Fb>::get_x()
{
    return 0;
//    return CGAL::circumcenter( this->vertex(0).point(),this->vertex(1).point(),this->vertex(2).point()).x();
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_y()
{
   return CGAL::circumcenter( this->vertex(0)->point(),this->vertex(1)->point(),this->vertex(2)->point()).y();
}

template < class Gt, class Fb>
double face<Gt, Fb>::get_z()
{
    return 0;
//    return CGAL::circumcenter( this->vertex(0).point(),this->vertex(1).point(),this->vertex(2).point()).y();
}

