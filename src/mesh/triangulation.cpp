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

triangulation::triangulation()
{
    LOG_WARNING << "No Matlab engine, plotting and all Matlab functionality will be disabled";
    #ifdef MATLAB
    _engine = NULL;
    _gfx = NULL;
#endif
    _size = 0;
}

#ifdef MATLAB
triangulation::triangulation(boost::shared_ptr<maw::matlab_engine> engine)
{
    _engine = engine;
    _gfx = boost::make_shared<maw::graphics>(_engine.get());
    _size = 0;
}
#endif
triangulation::~triangulation()
{
#ifdef MATLAB
    _engine = NULL;
    _gfx = NULL;
#endif
}

void triangulation::init(vector x, vector y, vector z)
{
    //    _mesh = boost::make_shared<Delaunay>();
    for (size_t i = 0; i < x->size(); i++)
    {
        Point p(x->at(i), y->at(i), z->at(i));
        Delaunay::Vertex_handle Vh = this->insert(p);
        Vh->set_id(i);
        ++i;
    }

    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(this->size()) + " triangles";
}

size_t triangulation::size()
{
    return _size;
}

mesh_elem triangulation::locate_face(double x, double y)
{
    Point_3 p(x, y, 0);

    Delaunay::Locate_type lt0;
    int li;

    mesh_elem m = this->locate(p, lt0, li);


    if (lt0 != Delaunay::FACE)
        return NULL;
    else
        return m;
}

#ifdef NOMATLAB
void triangulation::plot_time_series(double x, double y, std::string ID)
{
    if (!_engine)
    {
        LOG_WARNING << "No Matlab engine, plotting is disabled";
        return;
    }
    mesh_elem m = this->locate_face(x, y);

    if (m == NULL)
        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Couldn't find triangle at (x,y)"));

    maw::d_vec v(new arma::vec(m->face_time_series(ID).size()));

    time_series::variable_vec ts = m->face_time_series(ID);

    for (size_t i = 0; i < v->size(); i++)
    {
        (*v)(i) = ts.at(i);
    }

    _engine->put_double_matrix(ID, v);
    double handle = _gfx->plot_line(ID);
    _gfx->add_title(ID);
    _gfx->spin_until_close(handle);
}
#endif
void triangulation::from_file(std::string file)
{
    std::ifstream in(file);
    if (in.fail())
    {
        BOOST_THROW_EXCEPTION(file_read_error()
                << boost::errinfo_errno(errno)
                << boost::errinfo_file_name(file)
                );

    }

    Point pt;
    size_t i = 0;
    std::vector< K::Point_2 > pts;
    while (in >> pt)
    {
        Vertex_handle Vh = this->insert(pt);
        Vh->set_id(i);
        ++i;
        pts.push_back(K::Point_2(pt.x(), pt.y()));
    }
    _size = this->number_of_faces();
    _data_size = i;
    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(this->size()) + " triangles";

    _bbox = CGAL::bounding_box(pts.begin(), pts.end());
    //    std::cout << "0:"<< _bbox[0] << "1:"<< _bbox[1]<<"2:"<< _bbox[2]<<"3:"<< _bbox[3]<<std::endl;
    //    std::cout << _bbox[0] - _bbox[1] << std::endl;
    //    std::cout << _bbox[2] - _bbox[3] << std::endl;
}

void triangulation::to_file(double x, double y, std::string fname)
{
    mesh_elem m = this->locate_face(x, y);
    to_file(m, fname);
}

void triangulation::to_file(mesh_elem m, std::string fname)
{
    if (m == NULL)
        BOOST_THROW_EXCEPTION(mesh_error() << errstr_info("Couldn't find triangle at (x,y)"));

    m->to_file(fname);
}
#ifdef NOMATLAB

void triangulation::plot(std::string ID)
{
    if (!_engine)
    {
        LOG_WARNING << "No Matlab engine, plotting is disabled";
        return;
    }
    LOG_DEBUG << "Sending triangulation to matlab...";



    maw::d_mat tri(new arma::mat(_size, 3));
    maw::d_mat xyz(new arma::mat(_data_size, 3));
    maw::d_mat cdata(new arma::mat(_size, 1));

    size_t i = 0;

    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
            fit != this->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;

        (*tri)(i, 0) = face->vertex(0)->get_id() + 1; //+1 because matlab indexing starts at 1
        (*tri)(i, 1) = face->vertex(1)->get_id() + 1;
        (*tri)(i, 2) = face->vertex(2)->get_id() + 1;

        Delaunay::Triangle t = this->triangle(face);

        //        std::cout  << "xyz1:" << (*tri)(i, 0)-1 <<std::endl;
        (*xyz)((*tri)(i, 0) - 1, 0) = t[0].x();
        (*xyz)((*tri)(i, 0) - 1, 1) = t[0].y();
        (*xyz)((*tri)(i, 0) - 1, 2) = t[0].z();


        //        std::cout  << "xyz2:"<<(*tri)(i, 1)-1 <<std::endl;
        (*xyz)((*tri)(i, 1) - 1, 0) = t[1].x();
        (*xyz)((*tri)(i, 1) - 1, 1) = t[1].y();
        (*xyz)((*tri)(i, 1) - 1, 2) = t[1].z();

        //        std::cout  << "xyz3:"<<(*tri)(i, 2)-1 <<std::endl;
        (*xyz)((*tri)(i, 2) - 1, 0) = t[2].x();
        (*xyz)((*tri)(i, 2) - 1, 1) = t[2].y();
        (*xyz)((*tri)(i, 2) - 1, 2) = t[2].z();


        double d = fit->face_data(ID);
        (*cdata)(i) = d;
        //        std::cout << i <<std::endl; 
        ++i;
    }


    LOG_DEBUG << "Sending data to matlab...";
    _engine->put_double_matrix("tri", tri);
    _engine->put_double_matrix("elevation_data", xyz);
    _engine->put_double_matrix("face_data", cdata);

    double handle = _gfx->plot_patch("[elevation_data(:,1) elevation_data(:,2) elevation_data(:,3)]", "tri", "face_data(:)");
    _gfx->add_title(ID);

    _gfx->spin_until_close(handle);
    //    _engine->evaluate("save lol.mat");
    _engine->evaluate("clear tri elevation_data face_data");
}
#endif

void triangulation::to_vtu(std::string file_name)
{
    size_t i = 0;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkCellArray> triangles =
            vtkSmartPointer<vtkCellArray>::New();


    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
            fit != this->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;
        Delaunay::Triangle t = this->triangle(face);

        //        std::cout  << "xyz1:" << (*tri)(i, 0)-1 <<std::endl;
        points->InsertNextPoint(t[0].x(), t[0].y(), t[0].z());
        points->InsertNextPoint(t[1].x(), t[1].y(), t[1].z());
        points->InsertNextPoint(t[2].x(), t[2].y(), t[2].z());
    }

    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
            fit != this->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;

        vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();

        tri->GetPointIds()->SetId(face->vertex(0)->get_id(), 
                                    face->vertex(1)->get_id());

        triangles->InsertNextCell(tri);


        //        double d = fit->face_data(ID);
        //        (*cdata)(i) = d;
        //        //        std::cout << i <<std::endl; 
        //        ++i;
    }

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TRIANGLE, triangles);
    
    
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file_name.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();

}