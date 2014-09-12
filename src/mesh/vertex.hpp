#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

#include <boost/tuple/tuple.hpp>

/**
* \class vertex_info
* Arbitrary data may be added to the vertex via this class
*/
struct vertex_info
{
    virtual ~vertex_info(){};
};

/**
* \class ex_vertex
*
* The base CGAL vertex class is extended by ex_vertex to allow for add arbitrary vertex information as well as for incorperating the conecept of an ID. This ID is a global
* reference that allows for creating a Matlab or VTK based triangulation datastructure.
*/
template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class ex_vertex : public Vb
{
    typedef Vb Base;
public:
    vertex_info* info;
    typedef typename Vb::Vertex_handle Vertex_handle;
    typedef typename Vb::Face_handle Face_handle;
    typedef typename Vb::Point Point;

    template < typename TDS2 >
    struct Rebind_TDS
    {
        typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
        typedef ex_vertex<Gt, Vb2> Other;
    };
private:
    size_t _id;
public:
    ex_vertex();

    ex_vertex(const Point & p);

    ex_vertex(const Point & p, Face_handle f);

    ex_vertex(Face_handle f);

    /**
    * Sets the vertex to have a given id. Generally this would be the current x,y,z point read in from a file, for example.
    * \param id Vertex global id
    */
    void set_id(size_t id);

    /**
    * Returns the global id
    * \return global id
    */
    size_t get_id();

};

template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex() : Base()
{
    info = NULL;
}
template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex(const Point & p) : Base(p)
{
    info = NULL;
}
template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex(const Point & p, Face_handle f) : Base(f, p)
{
    info = NULL;
}
template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex(Face_handle f) : Base(f)
{
    info = NULL;
}

template < class Gt, class Vb>
void ex_vertex<Gt, Vb>::set_id(size_t id)
{
    _id = id;
}

template < class Gt, class Vb>
size_t ex_vertex<Gt, Vb>::get_id()
{
    return _id;
}

