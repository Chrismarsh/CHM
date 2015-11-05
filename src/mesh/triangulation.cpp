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


#include "triangulation.hpp"

triangulation::triangulation()
{
    LOG_WARNING << "No Matlab engine, plotting and all Matlab functionality will be disabled";
#ifdef MATLAB
    _engine = NULL;
    _gfx = NULL;
#endif
    _num_faces = 0;
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
    
    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
            fit != this->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;
        _faces.push_back(face);
    }
    
    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(this->size_faces()) + " triangles";
}

size_t triangulation::size_faces()
{
    return _num_faces;
}

size_t triangulation::size_vertex()
{
    return _num_vertex;
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

    timeseries::variable_vec ts = m->face_time_series(ID);

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
    std::vector< K::Point_2 > pts;  //TODO: clean this up
    while (in >> pt)
    {
        Vertex_handle Vh = this->insert(pt);
        Vh->set_id(i);
        ++i;
        pts.push_back(K::Point_2(pt.x(), pt.y()));
        _vertexes.push_back(Vh);
    }
    _num_faces = this->number_of_faces();
    _num_vertex = i;
    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(this->size_faces()) + " triangles";

    //_bbox = CGAL::bounding_box(pts.begin(), pts.end());
    //    std::cout << "0:"<< _bbox[0] << "1:"<< _bbox[1]<<"2:"<< _bbox[2]<<"3:"<< _bbox[3]<<std::endl;
    //    std::cout << _bbox[0] - _bbox[1] << std::endl;
    //    std::cout << _bbox[2] - _bbox[3] << std::endl;
    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
            fit != this->finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;
        _faces.push_back(face);
    }
}

Delaunay::Vertex_handle triangulation::vertex(size_t i)
{
    return _vertexes.at(i);
}

Delaunay::Face_handle triangulation::face(size_t i)
{
    return _faces.at(i);
}

void triangulation::timeseries_to_file(double x, double y, std::string fname)
{
    mesh_elem m = this->locate_face(x, y);
    timeseries_to_file(m, fname);
}

void triangulation::timeseries_to_file(mesh_elem m, std::string fname)
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

vtkSmartPointer<vtkUnstructuredGrid> triangulation::mesh_to_vtkUstructuredGrid()
{

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(this->_num_vertex);

    vtkSmartPointer<vtkCellArray> triangles =
            vtkSmartPointer<vtkCellArray>::New();

    typedef vtkFloatArray vtkFArray;

    std::map<std::string, vtkSmartPointer<vtkFArray> > data;

    //assume that all the faces have the same number of variables and the same types of variables
    //by this point this should be a fair assumption

    Delaunay::Finite_faces_iterator f = this->finite_faces_begin();
    auto variables = f->variables();
    for(auto& v: variables)
    {
        data[v] = vtkSmartPointer<vtkFArray>::New();
        data[v]->SetName(v.c_str());
    }

    //handle elevation/aspect/slope
    data["Elevation"] = vtkSmartPointer<vtkFArray>::New();
    data["Elevation"]->SetName("Elevation");

    data["Slope"] = vtkSmartPointer<vtkFArray>::New();
    data["Slope"]->SetName("Slope");

    data["Aspect"] = vtkSmartPointer<vtkFArray>::New();
    data["Aspect"]->SetName("Aspect");

//    auto test = vtkSmartPointer<vtkFArray>::New();
//    test->SetName("TEST");
//    test->SetNumberOfComponents(3);

//    data["Mean Curvature"] = vtkSmartPointer<vtkFloatArray>::New();
//    data["Mean Curvature"]->SetName("Mean Curvature");
//
//    data["Gauss Curvature"] = vtkSmartPointer<vtkFloatArray>::New();
//    data["Gauss Curvature"]->SetName("Gauss Curvature");

//    data["Liston Curvature"] = vtkSmartPointer<vtkFloatArray>::New();
//    data["Liston Curvature"]->SetName("Liston Curvature");

//    CGALMongeForm monge_form;
//    CGALMongeViaJet monge_fit;

//        int degree = 1;
//
//        std::vector<Point_3> pts;
//        pts.assign(myPoints.begin(),myPoints.end());
//        //LOG_DEBUG << pts.size();
//        monge_form = monge_fit(pts.begin() , pts.end(), 2,2);
//
//        double k1 = monge_form.principal_curvatures ( 0 );
//        double k2 = monge_form.principal_curvatures ( 1 );
//        curve =  0.5*(k1+k2);
//        data["Mean Curvature"]->InsertNextTuple1(curve);
////        curve =  (k1*k2);
////        data["Gauss Curvature"]->InsertNextTuple1(curve);
//
    int ii=0;
    for (Delaunay::Finite_faces_iterator fit = finite_faces_begin();
         fit != finite_faces_end(); ++fit)
    {
        Delaunay::Face_handle face = fit;
        Delaunay::Triangle t = this->triangle(face);

        vtkSmartPointer<vtkTriangle> tri =
                vtkSmartPointer<vtkTriangle>::New();

        tri->GetPointIds()->SetId(0, fit->vertex(0)->get_id());
        tri->GetPointIds()->SetId(1, fit->vertex(1)->get_id());
        tri->GetPointIds()->SetId(2, fit->vertex(2)->get_id());


        points->SetPoint(fit->vertex(0)->get_id(), t[0].x(), t[0].y(), t[0].z());
        points->SetPoint(fit->vertex(1)->get_id(), t[1].x(), t[1].y(), t[1].z());
        points->SetPoint(fit->vertex(2)->get_id(), t[2].x(), t[2].y(), t[1].z());

        triangles->InsertNextCell(tri);


        for(auto& v: variables)
        {
            double d = fit->face_data(v);
            data[v]->InsertNextTuple1(d);
        }
        data["Elevation"]->InsertNextTuple1(fit->get_z());
        data["Slope"]->InsertNextTuple1(fit->slope());
        data["Aspect"]->InsertNextTuple1(fit->aspect());

        //TODO: remove this hard coded test & add ability to do this for any
//        auto dir =*data["VW_dir"]->GetTuple(ii);
//        auto mag = *data["VW"]->GetTuple(ii);
//        test->InsertNextTuple3( -mag * sin(dir * 3.14159/180.0), -mag*cos(dir* 3.14159/180.0),fit->slope() * cos(fit->slope()*3.14159/180.0));
//        ii++;
    }



    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->SetPoints(points);
    unstructuredGrid->SetCells(VTK_TRIANGLE, triangles);
//    unstructuredGrid->GetPointData()->SetVectors(test);

    vtkSmartPointer<vtkGeometryFilter> geometryFilter =  vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(unstructuredGrid);
    geometryFilter->Update();
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>(geometryFilter->GetOutput());

//    vtkSmartPointer<vtkCurvatures> curvaturesFilter = vtkSmartPointer<vtkCurvatures>::New();
//
//    curvaturesFilter->SetInputData(polydata);
//    curvaturesFilter->SetCurvatureTypeToMean();
//    curvaturesFilter->Update();
//
//    unstructuredGrid->GetPointData()->AddArray(curvaturesFilter->GetOutput()->GetPointData()->GetScalars());


    for(auto& m : data)
    {
        unstructuredGrid->GetCellData()->AddArray(m.second);
    }

    return unstructuredGrid;
}
void triangulation::mesh_to_vtu(std::string file_name)
{

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = this->mesh_to_vtkUstructuredGrid();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file_name.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(unstructuredGrid);
#else
    writer->SetInputData(unstructuredGrid);
#endif
    writer->Write();

}

void triangulation::mesh_to_ascii(std::string file_name)
{
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = this->mesh_to_vtkUstructuredGrid();

    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(unstructuredGrid);
    geometryFilter->Update();
    geometryFilter->GetOutput()->GetPointData()->SetScalars(unstructuredGrid->GetCellData()->GetScalars());


    vtkSmartPointer<vtkTransform> flattener = vtkSmartPointer<vtkTransform>::New();
    flattener->Scale(1.0,1.0,0.0);
    flattener->Update();

    vtkSmartPointer<vtkTransformFilter> filt = vtkSmartPointer<vtkTransformFilter>::New();
    filt->SetInputData(geometryFilter->GetOutput());
    filt->SetTransform(flattener);
    filt->Update();

    double* bounds = unstructuredGrid->GetBounds();


    // Create a grid of points to interpolate over
    vtkSmartPointer<vtkPlaneSource> gridPoints = vtkSmartPointer<vtkPlaneSource>::New();

    size_t dx = 25; //(meters)
    size_t dy = 25;

    double distx = bounds[1] - bounds[0];
    double disty = bounds[3] - bounds[2];

    size_t gridSizeX = ceil(distx/dx);
    size_t gridSizeY = ceil(disty/dy);

    gridPoints->SetResolution(gridSizeX, gridSizeY); //number of cells, NOT! cell size.
    gridPoints->SetOrigin(bounds[0],  bounds[2], 0);
    gridPoints->SetPoint1(bounds[1],  bounds[2], 0);
    gridPoints->SetPoint2(bounds[0], bounds[3], 0);
    gridPoints->Update();


    // Perform the interpolation
    vtkSmartPointer<vtkProbeFilter> probeFilter =  vtkSmartPointer<vtkProbeFilter>::New();
    probeFilter->SetSourceData(filt->GetOutput());
    probeFilter->SetInputData(gridPoints->GetOutput());
    probeFilter->Update();

//    probeFilter->GetOutput()->GetCellData()->GetScalars(<#(const char*)name#>)

    vtkSmartPointer<vtkXMLPolyDataWriter> gridWriter =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    gridWriter->SetFileName ( file_name.c_str());
    gridWriter->SetInputData(probeFilter->GetOutput());
    gridWriter->Write();

}

 boost::shared_ptr<segmented_AABB> triangulation::AABB(size_t rows, size_t cols)
 {
    boost::shared_ptr<segmented_AABB> AABB = boost::make_shared<segmented_AABB>();
    AABB->make(this,rows,cols);
    return AABB;
     
 }


void segmented_AABB::make( triangulation* domain, size_t rows, size_t cols)
{
    int n_segments = cols;
    int n_v_segments = rows;

    arma::vec midpoint_bottom_x(n_segments);
    arma::vec midpoint_top_x(n_segments);

    arma::vec midpoint_bottom_y(n_segments);
    arma::vec midpoint_top_y(n_segments);


    this->n_rows = rows;
    this->n_cols = cols;

    std::vector<K::Point_2> bboxpt;
    bboxpt.reserve( domain->number_of_vertices());
    for(auto v = domain->finite_vertices_begin(); v!=domain->finite_vertices_end();v++ )
    {
        bboxpt.push_back(K::Point_2(v->point().x(),v->point().y()));
    }
        
    K::Iso_rectangle_2 rot_bbox = CGAL::bounding_box(bboxpt.begin(), bboxpt.end());

    //need to construct new axis aligned BBR
//    arma::u32 index;

    //left most pt
    double x_left = rot_bbox.vertex(0).x();
//    double y_left = rot_bbox.vertex(0).y();

    //right most pt
    double x_right = rot_bbox.vertex(2).x();
//    double y_right = rot_bbox.vertex(2).y();

    //bottom most pt
    double y_bottom = rot_bbox.vertex(1).y();
//    double x_bottom = rot_bbox.vertex(1).x();

    //top most pt
    double y_top =rot_bbox.vertex(3).y();
//    double x_top = rot_bbox.vertex(3).x();


    //horizontal step size
    double h_dx = (x_right - x_left) / n_segments;
    double v_dy = (y_top - y_bottom) / n_v_segments;


    //resize the row dimension
    m_grid.resize(n_v_segments);

    //start top left
    //loop over rows
    for (int i = 0; i < n_v_segments; i++)
    {
        //set the y coordinate to the current row top left coordinate
        double h_y = y_top - v_dy*i;
        double h_x = x_left;

        //set number of cols
        m_grid[i].resize(n_segments);


        //loop over columns
        for (int j = 0; j < n_segments; j++)
        {

            arma::mat* t = new arma::mat(5, 2);


            *t << h_x << h_y - v_dy << arma::endr // bottom left
                    << h_x + h_dx << h_y - v_dy << arma::endr //bottom right
                    << h_x + h_dx << h_y << arma::endr // top right
                    << h_x << h_y << arma::endr //top left
                    << h_x << h_y - v_dy << arma::endr; // bottom left


            m_grid[i][j] = new rect(t);
//            m_grid[i][j]->triangles.reserve( ceil(domain->number_of_vertices()/(rows*cols)) ); //estimate an approx number of triangles/rect. Try to cut down on repeat memory alloc.
            h_x = h_x + h_dx;
        }
    }

    
}

rect* segmented_AABB::get_rect(size_t row, size_t col)
{
    return m_grid[row][col];
}

segmented_AABB::segmented_AABB()
{

}

bool segmented_AABB::pt_in_rect(Delaunay::Vertex_handle v, rect* r)
{
    return pt_in_rect(v->point().x(), v->point().y(), r);
}

bool segmented_AABB::pt_in_rect(double x, double y, rect* r)
{
    for (int i = 0; i < 4; i++)
    {
        double x0 = r->coord->operator()(i, 0);
        double y0 = r->coord->operator()(i, 1);
        double x1 = r->coord->operator()(i + 1, 0);
        double y1 = r->coord->operator()(i + 1, 1);

        double pt = ((y - y0)*(x1 - x0) - (x - x0)*(y1 - y0));
        if (pt < 0)
            return false;
        else if (pt == 0)
            return false;
    }
    return true;
}

segmented_AABB::~segmented_AABB()
{
    for (auto it = m_grid.begin(); it != m_grid.end(); it++)
    {
        for (auto jt = it->begin(); jt != it->end(); jt++)
        {
            delete *jt;
        }
    }

}
