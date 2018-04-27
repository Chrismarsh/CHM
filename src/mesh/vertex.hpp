//
// Canadian Hydrological Model - The Canadian Hydrological Model (CHM) is a novel
// modular unstructured mesh based approach for hydrological modelling
// Copyright (C) 2018 Christopher Marsh
//
// This file is part of Canadian Hydrological Model.
//
// Canadian Hydrological Model is free software: you can redistribute it and/or
// modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Canadian Hydrological Model is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Canadian Hydrological Model.  If not, see
// <http://www.gnu.org/licenses/>.
//

#pragma once

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>

#include <boost/tuple/tuple.hpp>

#include "exception.hpp"

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
   // vertex_info* info;
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
    std::map<std::string,vertex_info* > _vertex_module_data;
public:
    ex_vertex();

    ex_vertex(const Point & p);

    ex_vertex(const Point & p, Face_handle f);

    ex_vertex(Face_handle f);

    ~ex_vertex();

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

    template<typename T>
    T*get_module_data(std::string module);

    template<typename T>
    T*make_module_data(std::string module);

    void set_module_data(std::string module, vertex_info *vi);

};

template < class Gt, class Vb >
ex_vertex<Gt, Vb>::~ex_vertex()
{

};

template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex() : Base()
{

    _id = 0;
}

template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex(const Point & p) : Base(p)
{
    _id = 0;

}
template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex(const Point & p, Face_handle f) : Base(f, p)
{
    _id = 0;

}
template < class Gt, class Vb>
ex_vertex<Gt, Vb>::ex_vertex(Face_handle f) : Base(f)
{
    _id = 0;

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

template < class Gt, class Vb>
template<typename T>
T* ex_vertex<Gt, Vb>::make_module_data(std::string module)
{

    auto it = _vertex_module_data.find(module);

    //we don't already have this, make a new one.
    if(it == _vertex_module_data.end())
    {
        T* vi = new T;
        _vertex_module_data[module] = vi;
    }

    return get_module_data<T>(module);
}


template < class Gt, class Vb>
template<typename T>
T* ex_vertex<Gt, Vb>::get_module_data(std::string module)
{

    auto it = _vertex_module_data.find(module);

    //we don't have
    if(it == _vertex_module_data.end())
    {
        BOOST_THROW_EXCEPTION(module_data_error() << errstr_info ("No data for module " + module));
    }

    return dynamic_cast<T*>(it->second);

}

template < class Gt, class Vb>
void ex_vertex<Gt, Vb>::set_module_data(std::string module, vertex_info *vi)
{
    _vertex_module_data[module] = vi;
}