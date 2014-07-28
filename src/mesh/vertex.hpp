#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

#include <boost/tuple/tuple.hpp>

struct vertex_info
{
    virtual ~vertex_info(){};
};

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

    void set_id(size_t id);
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

