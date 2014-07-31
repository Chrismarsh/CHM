
#include "terrain_shadows.hpp"

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
//http://www.rosshemsley.co.uk/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/

static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
    if ((x != n) && (x % (n / 100 + 1) != 0)) return;

    float ratio = x / (float) n;
    int c = ratio * w;

    std::cout << std::setw(3) << (int) (ratio * 100) << "% [";
    for (int x = 0; x < c; x++) std::cout << "=";
    for (int x = c; x < w; x++) std::cout << " ";
    std::cout << "]\r" << std::flush;
}

terrain_shadow::terrain_shadow(std::string ID)
{
    _provides->push_back("shadowed");
    _provides->push_back("z_prime");
    this->ID = ID;
    _parallel_type = parallel::domain;
    LOG_DEBUG << "Successfully instantiated module " << this->ID;

}

void terrain_shadow::run(mesh domain, boost::shared_ptr<global> global_param)
{
    double A = global_param->solar_az();
    double E = global_param->solar_el();

    if (global_param->solar_el() < 5)
    {
        for (triangulation::Finite_faces_iterator fit = domain->finite_faces_begin(); fit != domain->finite_faces_end(); ++fit)
        {
            //interpolate the station data to the current element
            triangulation::Face_handle face = fit;
            face->set_face_data("z_prime", 0); //unshadowed
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

    std::cout << K << std::endl;
    std::vector<triangulation::Face_handle> faces;

    //compute the rotation of each vertex
    for (triangulation::Finite_faces_iterator fit = domain->finite_faces_begin();
            fit != domain->finite_faces_end(); ++fit)
    {
        //interpolate the station data to the current element
        triangulation::Face_handle face = fit;

        //for each vertex
        for (int i = 0; i < 3; i++)
        {
            arma::vec coord(3);
            triangulation::Point p;
            // Vertex 0
            coord(0) = face->vertex(i)->point().x();
            coord(1) = face->vertex(i)->point().y();
            coord(2) = face->vertex(i)->point().z();
            std::cout.precision(15);
            //            std::cout << "orig: " <<  std::fixed <<face->vertex(i)->point() <<std::endl;

            coord = K*coord;
            p = triangulation::Point(coord(0), coord(1), coord(2));
            //            std::cout << "prj: " << std::fixed << p <<std::endl;

            if (!face->vertex(i)->info)
                face->vertex(i)->info = new vertex_flag();
            vertex_flag* vf = reinterpret_cast<vertex_flag*> (face->vertex(i)->info);
            vf->prj_vertex = p;
            vf->org_vertex = face->vertex(i)->point();

        }
    }

    //modify the underlying triangulation to reflect the rotated vertices
    for (triangulation::Finite_faces_iterator fit = domain->finite_faces_begin();
            fit != domain->finite_faces_end(); ++fit)
    {
        triangulation::Face_handle face = fit;

        //the idea here is as follows: we are iterating over each face, however a vertex may belong to
        // >1 face, so if we blindly modify a vertex, we may doubly or triply rotate a vertex. What this logic does
        // is only rotate vertexes we haven't visited, done by checking the bool vitied flag. Then, if we have not visited a vertex
        // we push back the vertex index as well as the original x.y.z point (as we need to undo the rotation at the end)
        for (int i = 0; i < 3; i++)
        {
            if (!face->vertex(i)->info)
                BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Null vertex info"));

            vertex_flag* vf = reinterpret_cast<vertex_flag*> (face->vertex(i)->info);
            face->vertex(i)->set_point(vf->prj_vertex);
        }

        //init memory but do nothing with it here
        module_shadow_face_info* tv = new module_shadow_face_info;
        face->info = tv;

        //save the iterator for the next step (sort))
        faces.push_back(face);
    }

    //sort descending
    std::sort(faces.begin(), faces.end(),
            [](triangulation::Face_handle fa, triangulation::Face_handle fb)->bool
            {
                return fa->center().z() > fb->center().z();
            });

    //            for(auto itr: faces)
    //            {
    //                std::cout << itr->center().z() << std::endl;
    //            }


    for (size_t j = 0; j < faces.size(); j++)
    {
        loadbar(j, faces.size());

        triangulation::Face_handle face_j = faces.at(j);

        K::Triangle_2 tj(K::Point_2(face_j->vertex(0)->point().x(), face_j->vertex(0)->point().y()),
                K::Point_2(face_j->vertex(1)->point().x(), face_j->vertex(1)->point().y()),
                K::Point_2(face_j->vertex(2)->point().x(), face_j->vertex(2)->point().y()));
        
        CGAL::Bbox_2 bj(tj.bbox());
        //compare to other triangles
        for (size_t k = j + 1; k < faces.size(); k++)
        {
            triangulation::Face_handle face_k = faces.at(k);
            module_shadow_face_info* face_k_info = reinterpret_cast<module_shadow_face_info*> (face_k->info);

            if (face_k_info->shadow == 0 && face_j->get_z() > face_k->get_z()) //tj is above tk, and tk is shadded by tj?
            {
                K::Triangle_2 tk(K::Point_2(face_k->vertex(0)->point().x(), face_k->vertex(0)->point().y()),
                        K::Point_2(face_k->vertex(1)->point().x(), face_k->vertex(1)->point().y()),
                        K::Point_2(face_k->vertex(2)->point().x(), face_k->vertex(2)->point().y()));

                CGAL::Bbox_2 bk(tk.bbox());

                if (CGAL::do_overlap(bk, bj))
                {
                    bool collision = face_k->intersects(face_j);
                    if (collision)
                    {
                        face_k_info->shadow = 1;
                    }
                }
            }
        }
    }

    for (triangulation::Finite_faces_iterator fit = domain->finite_faces_begin();
            fit != domain->finite_faces_end(); ++fit)
    {
        triangulation::Face_handle face = fit;

        face->set_face_data("z_prime", face->center().z());

        module_shadow_face_info* face_info = reinterpret_cast<module_shadow_face_info*> (face->info);
        face->set_face_data("shadowed", face_info->shadow);
    }

    // here we need to 'undo' the rotation we applied. It is a trade off of memor v. CPU. Frequently memory is the limitation
    // so we will do more work inorder to reuse our triangulation 
    for (triangulation::Finite_faces_iterator fit = domain->finite_faces_begin();
            fit != domain->finite_faces_end(); ++fit)
    {
        triangulation::Face_handle face = fit;
        for (int i = 0; i < 3; i++)
        {
            if (!face->vertex(i)->info)
                BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Null vertex info"));

            vertex_flag* vf = reinterpret_cast<vertex_flag*> (face->vertex(i)->info);
            face->vertex(i)->set_point(vf->org_vertex);
        }

    }
}

terrain_shadow::~terrain_shadow()
{


}



