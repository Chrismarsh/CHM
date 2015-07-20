
#include "Marsh_shading_iswr.hpp"

Marsh_shading_iswr::Marsh_shading_iswr(std::string ID)
{
    _provides->push_back("shadowed");
    _provides->push_back("z_prime");
    this->ID = ID;
    _parallel_type = parallel::domain;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;

}

void Marsh_shading_iswr::run(mesh domain, boost::shared_ptr<global> global_param)
{
    double A = global_param->solar_az();
    double E = global_param->solar_el();

    if (global_param->solar_el() < 5)
    {
#pragma omp parallel for
        for (size_t i = 0; i < domain->size_faces(); i++)
        {
            auto face = domain->face(i);
            face->set_face_data("z_prime", 0); //unshadowed
            face->set_face_data("shadowed", 0);
        }
        return;
    }
    //euler rotation matrix K
    arma::mat K;
    // eqns(6) & (7) in Montero
    double z0 = M_PI - A * M_PI / 180.0;
    double q0 = M_PI / 2.0 - E * M_PI / 180.0;

    K << cos(z0) << sin(z0) << 0 << arma::endr
            << -cos(q0) * sin(z0) << cos(q0) * cos(z0) << sin(q0) << arma::endr
            << sin(q0) * sin(z0) << -cos(z0) * sin(q0) << cos(q0) << arma::endr;


    //compute the rotation of each vertex

    //    tbb::concurrent_vector<triangulation::Face_handle> rot_faces;
    //    rot_faces.grow_by(domain->size());

#pragma omp parallel for
    for (size_t i = 0; i < domain->size_vertex(); i++)
    {
        auto vert = domain->vertex(i);
        if (!vert->info)
            vert->info = new vertex_data();
        vertex_data * vf = reinterpret_cast<vertex_data *> (vert->info);

        arma::vec coord(3);
        triangulation::Point p;

        coord(0) = vert->point().x();
        coord(1) = vert->point().y();
        coord(2) = vert->point().z();

        coord = K*coord;
        p = triangulation::Point(coord(0), coord(1), coord(2));
        vf->prj_vertex = p;
        vf->org_vertex = vert->point();

        vert->set_point(vf->prj_vertex); //modify the underlying triangulation to reflect the rotated vertices
    }

    auto BBR = domain->AABB(10,10);

    #pragma omp parallel for
    for (size_t i = 0; i < domain->size_faces(); i++)
    {
        auto face = domain->face(i);

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
        module_shadow_face_info* tv = new module_shadow_face_info;
//        face->info = tv;
        face->set_module_face_data(ID, tv);
        tv->z_prime = CGAL::centroid(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point()).z();
    }

    LOG_DEBUG << "AABB is " <<BBR->n_rows << "x" << BBR->n_rows;
    for (size_t j = 0; j < BBR->n_rows; j++)
    {
        for (size_t k = 0; k < BBR->n_cols; k++)
        {
            LOG_DEBUG << BBR->get_rect(j, k)->triangles.size();
        }
    }

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
//                        module_shadow_face_info* fa_info = reinterpret_cast<module_shadow_face_info*> (fa->info);
//                        module_shadow_face_info* fb_info = reinterpret_cast<module_shadow_face_info*> (fb->info);
                        module_shadow_face_info* fa_info = reinterpret_cast<module_shadow_face_info*> (fa->module_face_data("shadows"));
                        module_shadow_face_info* fb_info = reinterpret_cast<module_shadow_face_info*> (fb->module_face_data("shadows"));

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
                    module_shadow_face_info* face_k_info = reinterpret_cast<module_shadow_face_info*> (face_k->module_face_data(ID));

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



                module_shadow_face_info* face_info = reinterpret_cast<module_shadow_face_info*> (face_j->module_face_data(ID));
                face_j->set_face_data("shadowed", face_info->shadow);
                face_j->set_face_data("z_prime", face_info->z_prime);


            }
        }
    }

    // here we need to 'undo' the rotation we applied.
#pragma omp parallel for
    for (size_t i = 0; i < domain->size_vertex(); i++)
    {
        auto vert = domain->vertex(i);
        if (!vert->info)
            BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Null vertex info"));

        vertex_data * vf = reinterpret_cast<vertex_data *> (vert->info);
        vert->set_point(vf->org_vertex);


    }
}

Marsh_shading_iswr::~Marsh_shading_iswr()
{


}



