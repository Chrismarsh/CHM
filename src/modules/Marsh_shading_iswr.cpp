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


#include "Marsh_shading_iswr.hpp"
REGISTER_MODULE_CPP(Marsh_shading_iswr);

Marsh_shading_iswr::Marsh_shading_iswr(config_file cfg)
        :module_base("Marsh_shading_iswr", parallel::domain, cfg)
{
    depends("solar_az");
    depends("solar_el");

    provides("shadow");
    provides("z_prime");

    x_AABB = cfg.get<int>("x_AABB",10);
    y_AABB = cfg.get<int>("y_AABB",10);
    LOG_DEBUG << "Successfully instantiated module " << this->ID;

}

void Marsh_shading_iswr::run(mesh& domain)
{
    //compute the rotation of each vertex

    //    tbb::concurrent_vector<triangulation::Face_handle> rot_faces;
    //    rot_faces.grow_by(domain->size());

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_vertex(); i++)
    {
        auto vert = domain->vertex(i);

        double A =  vert->face()->face_data("solar_az");
        double E = vert->face()->face_data("solar_el");


        //euler rotation matrix K
        arma::mat K;
        // eqns(6) & (7) in Montero
        double z0 = M_PI - A * M_PI / 180.0;
        double q0 = M_PI / 2.0 - E * M_PI / 180.0;

        K << cos(z0) << sin(z0) << 0 << arma::endr
          << -cos(q0) * sin(z0) << cos(q0) * cos(z0) << sin(q0) << arma::endr
          << sin(q0) * sin(z0) << -cos(z0) * sin(q0) << cos(q0) << arma::endr;


        auto vf = vert->make_module_data<vertex_data>(ID);

        arma::vec coord(3);

        Point_3 p;
        coord(0) = vert->point().x();
        coord(1) = vert->point().y();
        coord(2) = vert->point().z();

        coord = K*coord;
        p = Point_3(coord(0), coord(1), coord(2));
        vf->prj_vertex = p;
        vf->org_vertex = vert->point();

        vert->set_point(vf->prj_vertex); //modify the underlying triangulation to reflect the rotated vertices

    }

    auto BBR = domain->AABB(x_AABB,y_AABB);

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);
        double E = domain->face(i)->face_data("solar_el");
        if (E < 5)
        {
            (*face)["z_prime"_s]= 0; //unshadowed
            (*face)["shadow"_s]= 0;
            continue;
        }

        for (size_t j = 0; j < BBR->n_rows; j++)
        {
            for (size_t k = 0; k < BBR->n_cols; k++)
            {

                if (BBR->pt_in_rect(face->vertex(0), BBR->get_rect(j, k)) || //pt1
                        BBR->pt_in_rect(face->vertex(1), BBR->get_rect(j, k)) || //pt2
                        BBR->pt_in_rect(face->vertex(2), BBR->get_rect(j, k))) //pt3
                {
//#pragma omp critical
//                    {
                        BBR->get_rect(j, k)->triangles.push_back(face);
//                    }
                }
            }

        }
        //init memory
        auto tv = face->make_module_data<module_shadow_face_info>(ID);
       // module_shadow_face_info* tv = new module_shadow_face_info;
        //face->set_module_data(ID, tv);
        tv->z_prime = CGAL::centroid(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point()).z();
    }

//    LOG_DEBUG << "AABB is " <<BBR->n_rows << "x" << BBR->n_rows;
//    for (size_t j = 0; j < BBR->n_rows; j++)
//    {
//        for (size_t k = 0; k < BBR->n_cols; k++)
//        {
//            LOG_DEBUG << BBR->get_rect(j, k)->triangles.size();
//        }
//    }

    tbb::task_scheduler_init init;

#pragma omp parallel for
    for (size_t i = 0; i < BBR->n_rows; i++)
    {
        for (size_t ii = 0; ii < BBR->n_cols; ii++)
        {
            //sort descending
            tbb::parallel_sort(BBR->get_rect(i, ii)->triangles.begin(), BBR->get_rect(i, ii)->triangles.end(),
                    [](triangulation::Face_handle fa, triangulation::Face_handle fb)->bool
                    {
                          auto fa_info = fa->get_module_data<module_shadow_face_info>(
                                  "Marsh_shading_iswr");
                          auto fb_info = fb->get_module_data<module_shadow_face_info>(
                                  "Marsh_shading_iswr");

                        return fa_info->z_prime > fb_info->z_prime;
                    });

            size_t num_tri = BBR->get_rect(i, ii)->triangles.size();


            for (size_t j = 0; j < num_tri; j++)
            {
                triangulation::Face_handle face_j = BBR->get_rect(i, ii)->triangles.at(j);

                K::Triangle_2 tj(K::Point_2(face_j->vertex(0)->point().x(), face_j->vertex(0)->point().y()),
                        K::Point_2(face_j->vertex(1)->point().x(), face_j->vertex(1)->point().y()),
                        K::Point_2(face_j->vertex(2)->point().x(), face_j->vertex(2)->point().y()));

                CGAL::Bbox_2 bj(tj.bbox());
                //compare to other triangles

                for (size_t k = j + 1; k < num_tri; k++)
                {
                    triangulation::Face_handle face_k = BBR->get_rect(i, ii)->triangles.at(k);
                     auto face_k_info = face_k->get_module_data<module_shadow_face_info>(ID);

                    if (face_k_info->shadow == 0)
                    {
                        //not needed as our sort will ensure face_j > face_k
                        if(face_j->get_z() > face_k->get_z())  //tj is above tk, and tk is shadded by tj?)
                        {
                            K::Triangle_2 tk(K::Point_2(face_k->vertex(0)->point().x(), face_k->vertex(0)->point().y()),
                                    K::Point_2(face_k->vertex(1)->point().x(), face_k->vertex(1)->point().y()),
                                    K::Point_2(face_k->vertex(2)->point().x(), face_k->vertex(2)->point().y()));

                            CGAL::Bbox_2 bk(tk.bbox());

                            if (CGAL::do_overlap(bk, bj)) {
                                bool collision = face_k->intersects(face_j);
                                if (collision) {

                                    face_k_info->shadow = 1;
                                }
                            }
                        }
                    }
                }

               auto face_info = face_j->get_module_data<module_shadow_face_info>(ID);
                face_j->set_face_data("shadow", face_info->shadow);
                face_j->set_face_data("z_prime", face_info->z_prime);

            }
        }
    }

    // here we need to 'undo' the rotation we applied.
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_vertex(); i++)
    {
        auto vert = domain->vertex(i);
        auto  vf = vert->get_module_data<vertex_data>(ID);
        vert->set_point(vf->org_vertex);
    }

}

Marsh_shading_iswr::~Marsh_shading_iswr()
{


}
