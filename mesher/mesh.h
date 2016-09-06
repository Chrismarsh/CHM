#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include "mesh_2_criteria_area.h"

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

    typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K> Vb;

    template<class Gt, class Fb = CGAL::Delaunay_mesh_face_base_2<Gt> >
    class Delaunay_mesh_face_base_info_2 : public Fb
    {
    public:
        typedef Gt Geom_traits;
        typedef typename Fb::Vertex_handle Vertex_handle;
        typedef typename Fb::Face_handle Face_handle;

        int id;

        template<typename TDS2>
        struct Rebind_TDS
        {
            typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
            typedef Delaunay_mesh_face_base_info_2<Gt, Fb2> Other;
        };

        Delaunay_mesh_face_base_info_2() : Fb()
        {
            id = 0;
        }

        Delaunay_mesh_face_base_info_2(Vertex_handle v0,
                                       Vertex_handle v1,
                                       Vertex_handle v2)
                : Fb(v0, v1, v2)
        {
            id = 0;
        }


        Delaunay_mesh_face_base_info_2(Vertex_handle v0,
                                       Vertex_handle v1,
                                       Vertex_handle v2,
                                       Face_handle n0,
                                       Face_handle n1,
                                       Face_handle n2)
                : Fb(v0, v1, v2, n0, n1, n2)
        {
            id = 0;
        }


    };

    typedef Delaunay_mesh_face_base_info_2<K> Fb;

    typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
    typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

    typedef mesh_2_criterion_area<CDT> Criteria;
    typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

    typedef CDT::Vertex_handle Vertex_handle;
    typedef CDT::Face_handle Face_handle;
    typedef CDT::Face_circulator Face_circulator;
    typedef CDT::Point Point;
