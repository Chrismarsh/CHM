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


#include "triangulation.h"

triangulation::triangulation(boost::shared_ptr<maw::matlab_engine> engine)
{
    _engine = engine;
    _gfx = boost::make_shared<maw::graphics>(_engine.get());
}

triangulation::~triangulation()
{
    _engine = NULL;
    _gfx = NULL;
}

void triangulation::init(vector x, vector y, vector z)
{
    if (!_engine)
    {
        BOOST_THROW_EXCEPTION(mesh_error()
                << errstr_info("No matlab engine"));
    }

    //    _mesh = boost::make_shared<Delaunay>();
    for (size_t i = 0; i < x->size(); i++)
    {
        Point p(x->at(i), y->at(i), z->at(i));
        Delaunay::Vertex_handle Vh = this->insert(p);
        Vh->set_id(i);
        ++i;
    }

    //    _mesh
    //    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(m_size) +" triangles";
}

size_t triangulation::size()
{
    return _size;
}



//mesh_elem triangulation::operator()(size_t t)
//{
//    Delaunay::Finite_faces_iterator fit = _mesh->finite_faces_begin();
//    for(size_t i=0;i<t;i++)
//        fit++;
//    
//    Delaunay::Face_handle face = fit;
//    return *face;
//}

void triangulation::set_vertex_data(vector data)
{


}

//triangle* triangulation::find_containing_triangle(double x, double y)
//{
////    for (std::vector<triangle*>::iterator it = m_triangles.begin(); it != m_triangles.end(); it++)
////    {
////        if ((*it)->contains(x, y))
////            return *it;
////    }
//
//    return NULL;
//}



//maw::d_vec triangulation::mf_face_data(std::string ID)
//{
//    maw::d_vec data(new arma::vec(m_size));
//
//    int i = 0;
//    for (auto& it : m_triangles)
//    {
//        double d = it->get_face_data(ID);
//        (*data)(i) = d;
//        ++i;
//
//    }
//
//    return data;
//}

//void triangulation::plot_time_series(double x, double y, std::string ID)
//{
//    mesh_elem* m  = _mesh->find_containing_triangle(x,y);
//    
//    if(!m)
//        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Couldn't find triangle at (x,y)"));
//    
//    maw::d_vec v(new arma::vec(m->get_face_time_series(ID).size()));
//    
////    v->resize(m->get_face_time_series(ID).size());
//    
//    for(size_t i=0; i < v->size(); i++ )
//    {
//        (*v)(i) = m->get_face_time_series(ID).at(i);
//    }
//    
//    _engine->put_double_matrix(ID,v);
//   double handle = _gfx->plot_line(ID);
//    _gfx->spin_until_close(handle);
//}

void triangulation::from_file(std::string file)
{
    std::ifstream in(file);
    Point pt;
    size_t i = 0;
    while (in >> pt)
    {
        Vertex_handle Vh = this->insert(pt);
        Vh->set_id(i);
        ++i;
    }
    _size = this->number_of_faces();

}

void triangulation::plot(std::string ID)
{
    LOG_DEBUG << "Sending triangulation to matlab...";

    maw::d_mat tri(new arma::mat(_size, 3));
    maw::d_mat xyz(new arma::mat(_size, 3));
    maw::d_mat face(new arma::mat(_size, 1));

    size_t i = 0;

    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
            fit != this->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;

        (*tri)(i, 0) = face->vertex(0)->get_id() + 1;
        (*tri)(i, 1) = face->vertex(1)->get_id() + 1;
        (*tri)(i, 2) = face->vertex(2)->get_id() + 1;

        Delaunay::Triangle t = this->triangle(face);

        (*xyz)((*tri)(i, 0), 0) = t[0].x();
        (*xyz)((*tri)(i, 0), 1) = t[0].y();
        (*xyz)((*tri)(i, 0), 2) = t[0].z();

        (*xyz)((*tri)(i, 1), 0) = t[1].x();
        (*xyz)((*tri)(i, 1), 1) = t[1].y();
        (*xyz)((*tri)(i, 1), 2) = t[1].z();

        (*xyz)((*tri)(i, 2), 0) = t[2].x();
        (*xyz)((*tri)(i, 2), 1) = t[2].y();
        (*xyz)((*tri)(i, 2), 2) = t[2].z();

        //        double d = it->get_face_data(ID);
        //        (*data)(i) = d;
        //        ++i;
    }


    LOG_DEBUG << "Sending data to matlab...";
    _engine->put_double_matrix("tri", tri);
    _engine->put_double_matrix("elevation_data", xyz);
    _engine->put_double_matrix("face_data", face);

    double handle = _gfx->plot_patch("[elevation_data(:,1) elevation_data(:,2) elevation_data(:,3)]", "tri", "face_data(:)");
    _gfx->add_title(ID);

    _gfx->spin_until_close(handle);
    //    _engine->evaluate("save lol.mat");
    _engine->evaluate("clear tri elevation_data face_data");
}


