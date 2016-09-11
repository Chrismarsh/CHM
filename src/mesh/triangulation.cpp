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
 //   LOG_WARNING << "No Matlab engine, plotting and all Matlab functionality will be disabled";
#ifdef MATLAB
    _engine = NULL;
    _gfx = NULL;
#endif
    _num_faces = 0;
    _vtk_unstructuredGrid = nullptr;
    _is_geographic = false;
    _UTM_zone = 0;
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

//void triangulation::init(vector x, vector y, vector z)
//{
//    //    _mesh = boost::make_shared<Delaunay>();
//    for (size_t i = 0; i < x->size(); i++)
//    {
//        Point_3 p(x->at(i), y->at(i), z->at(i));
//        Delaunay::Vertex_handle Vh = this->insert(p);
//        Vh->set_id(i);
//        ++i;
//    }
//
//    for (Delaunay::Finite_faces_iterator fit = this->finite_faces_begin();
//            fit != this->finite_faces_end(); ++fit)
//    {
//        Delaunay::Face_handle face = fit;
//        _faces.push_back(face);
//    }
//
//    LOG_DEBUG << "Created a mesh with " + boost::lexical_cast<std::string>(this->size_faces()) + " triangles";
//}

std::string triangulation::wkt()
{
    return _srs_wkt;
}
int triangulation::UTM_zone()
{
    return _UTM_zone;
}
bool triangulation::is_geographic()
{
    return _is_geographic;
}
size_t triangulation::size_faces()
{
    return _num_faces;
}

size_t triangulation::size_vertex()
{
    return _num_vertex;
}

mesh_elem triangulation::locate_face(Point_2 query)
{
    K_neighbor_search search(*(dD_tree.get()), query, 1);
    auto it = search.begin();
    return boost::get<1>(it->first);
}
mesh_elem triangulation::locate_face(double x, double y)
{
    //http://doc.cgal.org/latest/Spatial_searching/index.html
    Point_2 query(x,y);

    return locate_face(query);

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

void triangulation::serialize_parameter(std::string output_path, std::string parameter)
{

    pt::ptree out;
    pt::ptree subtree;

    for(auto& face : _faces)
    {
        pt::ptree item;
        item.put_value(face->get_parameter(parameter));
        subtree.push_back(std::make_pair("",item));
    }
    out.put_child( pt::ptree::path_type(parameter),subtree);

    pt::write_json(output_path,out);

}
std::set<std::string>  triangulation::from_json(pt::ptree &mesh)
{

    size_t nvertex_toread = mesh.get<size_t>("mesh.nvertex");
    LOG_DEBUG << "Reading in #vertex=" << nvertex_toread;
    size_t i=0;

    //paraview struggles with lat/long as it doesn't seem to have the accuracy. So we need to scale up lat-long.
    int is_geographic = mesh.get<int>("mesh.is_geographic");
    if( is_geographic == 1)
    {
        _is_geographic = true;
    }


    if(!_is_geographic)
    {
        _UTM_zone = mesh.get<int>("mesh.UTM_zone");
    }

//    _srs_wkt = mesh.get<std::string>("mesh.srs_wkt");

    for (auto &itr : mesh.get_child("mesh.vertex"))
    {
        std::vector<double> items;
        //iterate over the vertex triples
        for(auto& jtr: itr.second)
        {
            items.push_back(jtr.second.get_value<double>());
        }
        Point_3 pt( items[0], items[1], items[2]);

        Vertex_handle Vh = this->create_vertex();
        Vh->set_point(pt);
        Vh->set_id(i);
        _vertexes.push_back(Vh);
        i++;
    }
    _num_vertex = this->number_of_vertices();
    LOG_DEBUG << "# nodes created = " << _num_vertex;

    if( this->number_of_vertices() != nvertex_toread)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                "Expected: " + std::to_string(nvertex_toread) + " vertex, got: " + std::to_string(this->number_of_vertices())));
    }

    //read in faces
    size_t num_elem = mesh.get<int>("mesh.nelem");
    LOG_DEBUG << "Reading in #elem = " << num_elem;
    //set our mesh dimensions
    this->set_dimension(2);

    //vectors to hold the center of a face, and a pointer to that face to generate the spatial search tree
    std::vector<Point_2> center_points;

    i = 0;
    size_t cid = 0;
    for (auto &itr : mesh.get_child("mesh.elem"))
    {
        std::vector<int> items;
        //iterate over the vertex triples
        for(auto& jtr: itr.second)
        {
            items.push_back(jtr.second.get_value<size_t>());
        }
        auto vert1 = _vertexes.at(items[0]); //0 indexing
        auto vert2 = _vertexes.at(items[1]);
        auto vert3 = _vertexes.at(items[2]);

        auto face = this->create_face(vert1,vert2,vert3);
        face->cell_id = cid++;
        if( is_geographic == 1)
        {
            face->_is_geographic = true;
        }
        face->_debug_ID= --i; //all ids will be negative starting at -1. Named ids (for output) will be positive starting at 0
        face->_debug_name= std::to_string(i);
        _faces.push_back(face);

        Point_2 pt2(face->center().x(),face->center().y());

        center_points.push_back(pt2);
    }

    //make the search tree
    dD_tree = boost::make_shared<Tree>(boost::make_zip_iterator(boost::make_tuple( center_points.begin(),_faces.begin() )),
                                 boost::make_zip_iterator(boost::make_tuple( center_points.end(),  _faces.end() ) )
    );

    _num_faces = this->number_of_faces();

    LOG_DEBUG << "Created a mesh with " << this->size_faces() << " triangles";

    if( this->number_of_faces() != num_elem)
    {
        BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                "Expected: " + std::to_string(num_elem) + " elems, got: " + std::to_string(this->size_faces())));
    }

    LOG_DEBUG << "Building face neighbours";
    i=0;
    int nelem = mesh.get<int>("mesh.nelem"); // what we are expecting to see, 0 indexed
    for (auto &itr : mesh.get_child("mesh.neigh"))
    {
        std::vector<int> items;
        //iterate over the vertex triples
        for(auto& jtr: itr.second)
        {
            items.push_back(jtr.second.get_value<int>());
        }

        auto face = _faces.at(i);

        if(    items[0] > nelem
            || items[1] > nelem
            || items[2] > nelem)
        {
            BOOST_THROW_EXCEPTION(config_error() << errstr_info(
                    "Face " + std::to_string(i) + " has out of bound neighbours."));
        }
        //-1 is now the no neighbour value
        Face_handle face0 =  items[0] != -1 ?_faces.at( items[0] ) : Face_handle();
        Face_handle face1 =  items[1] != -1 ?_faces.at( items[1] ) : Face_handle();
        Face_handle face2 =  items[2] != -1 ?_faces.at( items[2] ) : Face_handle();

        face->set_neighbors(face0,face1,face2);

        i++;
    }
    //loop over parameter name

    std::set<std::string> parameters;
    try
    {
        for (auto &itr : mesh.get_child("parameters"))
        {
            i = 0; // reset evertime we get a new parameter set
            auto name = itr.first.data();
            LOG_DEBUG << "Applying parameter: " << name;
            for (auto &jtr : itr.second)
            {
                auto face = _faces.at(i);
                double value = jtr.second.get_value<double>();
                //            value == -9999. ? value = nan("") : value;
                face->set_parameter(name, value);
                i++;
            }

            parameters.insert(name);
        }
    }catch(pt::ptree_bad_path& e)
    {
        // we don't have this section, no worries
    }

    std::set<std::string> ics;
    try
    {
        for (auto &itr : mesh.get_child("initial_conditions"))
        {
            i=0; // reset evertime we get a new parameter set
            auto name = itr.first.data();
            LOG_DEBUG << "Applying IC: " << name;
            for (auto &jtr : itr.second)
            {
                auto face = _faces.at(i);
                double value = jtr.second.get_value<double>();
    //            alue == -9999. ? value = nan("") : value;
                face->set_initial_condition(name,value);
                i++;
            }

            ics.insert(name);
        }
    }catch(pt::ptree_bad_path& e)
    {
        // we don't have this section, no worries
    }
    _num_faces = this->number_of_faces();
    _num_vertex = this->number_of_vertices();

    return parameters;
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

void triangulation::init_vtkUnstructured_Grid(std::vector<std::string> output_variables)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(this->_num_vertex);
    //   points->Allocate(this->_num_vertex);

    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    triangles->Allocate(this->_num_vertex);

    double scale = is_geographic() == true ? 100000. : 1.;

    for (size_t i = 0; i < this->size_faces(); i++)
    {
        Delaunay::Face_handle fit = this->face(i);

        vtkSmartPointer<vtkTriangle> tri =
                vtkSmartPointer<vtkTriangle>::New();

        tri->GetPointIds()->SetId(0, fit->vertex(0)->get_id());
        tri->GetPointIds()->SetId(1, fit->vertex(1)->get_id());
        tri->GetPointIds()->SetId(2, fit->vertex(2)->get_id());

        points->SetPoint(fit->vertex(0)->get_id(), face(i)->vertex(0)->point().x()*scale, face(i)->vertex(0)->point().y()*scale, face(i)->vertex(0)->point().z());
        points->SetPoint(fit->vertex(1)->get_id(), face(i)->vertex(1)->point().x()*scale, face(i)->vertex(1)->point().y()*scale, face(i)->vertex(1)->point().z());
        points->SetPoint(fit->vertex(2)->get_id(), face(i)->vertex(2)->point().x()*scale, face(i)->vertex(2)->point().y()*scale, face(i)->vertex(2)->point().z());

        triangles->InsertNextCell(tri);
    }

    _vtk_unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    _vtk_unstructuredGrid->SetPoints(points);
    _vtk_unstructuredGrid->SetCells(VTK_TRIANGLE, triangles);



    //assume that all the faces have the same number of variables and the same types of variables
    //by this point this should be a fair assumption

//    Delaunay::Finite_faces_iterator f = this->face(0); //this->finite_faces_begin();


    auto variables = output_variables.size() == 0 ? this->face(0)->variables() : output_variables;
    for(auto& v: variables)
    {
        data[v] = vtkSmartPointer<vtkFloatArray>::New();
        data[v]->SetName(v.c_str());
    }

    auto params = this->face(0)->parameters();
    for(auto& v: params)
    {
        data["[param] " + v] = vtkSmartPointer<vtkFloatArray>::New();
        data["[param] " + v]->SetName(("[param] " + v).c_str());
    }

    auto ics = this->face(0)->initial_conditions();
    for(auto& v: ics)
    {
        data["[ic] " + v] = vtkSmartPointer<vtkFloatArray>::New();
        data["[ic] " + v]->SetName(("[ic] " + v).c_str());
    }


    //handle elevation/aspect/slope
    data["Elevation"] = vtkSmartPointer<vtkFloatArray>::New();
    data["Elevation"]->SetName("Elevation");

    data["Slope"] = vtkSmartPointer<vtkFloatArray>::New();
    data["Slope"]->SetName("Slope");

    data["Aspect"] = vtkSmartPointer<vtkFloatArray>::New();
    data["Aspect"]->SetName("Aspect");

    data["Area"] = vtkSmartPointer<vtkFloatArray>::New();
    data["Area"]->SetName("Area");

}

void triangulation::update_vtk_data(std::vector<std::string> output_variables)
{
    //if we haven't inited yet, do so.
    if(!_vtk_unstructuredGrid)
    {
        this->init_vtkUnstructured_Grid(output_variables);
    }



//    auto test = vtkSmartPointer<vtkFloatArray>::New();
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

    auto variables = output_variables.size() == 0 ? this->face(0)->variables() : output_variables;
    auto params = this->face(0)->parameters();
    auto ics = this->face(0)->initial_conditions();

    for (size_t i = 0; i < this->size_faces(); i++)
    {
        Delaunay::Face_handle fit = this->face(i);

        for (auto &v: variables)
        {
            double d = fit->face_data(v);
            if(d == -9999.)
            {
                d = nan("");
            }

            data[v]->InsertTuple1(i,d);
        }
        for(auto& v: params)
        {
            double d = fit->get_parameter(v);
            if(d == -9999.)
            {
                d = nan("");
            }
            data["[param] " +  v]->InsertTuple1(i,d);
        }

        for(auto& v: ics)
        {
            double d = fit->get_initial_condition(v);
            if(d == -9999.)
            {
                d = nan("");
            }
            data["[ic] "+ v]->InsertTuple1(i,d);
        }

        data["Elevation"]->InsertTuple1(i,fit->get_z());
        data["Slope"]->InsertTuple1(i,fit->slope());
        data["Aspect"]->InsertTuple1(i,fit->aspect());
        data["Area"]->InsertTuple1(i,fit->get_area());

//        data["Elevation"]->InsertNextTuple1(fit->get_z());
//        data["Slope"]->InsertNextTuple1(fit->slope());
//        data["Aspect"]->InsertNextTuple1(fit->aspect());
//        data["Area"]->InsertNextTuple1(fit->get_area());

        //TODO: remove this hard coded test & add ability to do this for any
//        auto dir =*data["VW_dir"]->GetTuple(ii);
//        auto mag = *data["VW"]->GetTuple(ii);
//        test->InsertNextTuple3( -mag * sin(dir * 3.14159/180.0), -mag*cos(dir* 3.14159/180.0),fit->slope() * cos(fit->slope()*3.14159/180.0));
//        ii++;
    }

    //    _vtk_unstructuredGrid->GetPointData()->SetVectors(test);

    for(auto& m : data)
    {
        _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);
    }

//    return _vtk_unstructuredGrid;
}
void triangulation::write_vtu(std::string file_name)
{
    //this now needs to be called from outside these functions
//    update_vtk_data();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file_name.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(_vtk_unstructuredGrid);
#else
    writer->SetInputData(_vtk_unstructuredGrid);
#endif
    writer->Write();

}

void triangulation::write_vtp(std::string file_name)
{
    //this now needs to be called from outside these functions
//    update_vtk_data();

//    vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
//    geometryFilter->SetInputData(_vtk_unstructuredGrid);
//    geometryFilter->Update();
//    geometryFilter->GetOutput()->GetPointData()->SetScalars(_vtk_unstructuredGrid->GetCellData()->GetScalars());
//
//
//    vtkSmartPointer<vtkTransform> flattener = vtkSmartPointer<vtkTransform>::New();
//    flattener->Scale(1.0,1.0,0.0);
//    flattener->Update();
//
//    vtkSmartPointer<vtkTransformFilter> filt = vtkSmartPointer<vtkTransformFilter>::New();
//    filt->SetInputData(geometryFilter->GetOutput());
//    filt->SetTransform(flattener);
//    filt->Update();
//
//    double* bounds = _vtk_unstructuredGrid->GetBounds();
//
//
//    // Create a grid of points to interpolate over
//    vtkSmartPointer<vtkPlaneSource> gridPoints = vtkSmartPointer<vtkPlaneSource>::New();
//
//    size_t dx = 300; //(meters)
//    size_t dy = 300;
//
//    double distx = bounds[1] - bounds[0];
//    double disty = bounds[3] - bounds[2];
//
//    size_t gridSizeX = ceil(distx/dx);
//    size_t gridSizeY = ceil(disty/dy);
//
//    gridPoints->SetResolution(gridSizeX, gridSizeY); //number of cells, NOT! cell size.
//    gridPoints->SetOrigin(bounds[0],  bounds[2], 0);
//    gridPoints->SetPoint1(bounds[1],  bounds[2], 0);
//    gridPoints->SetPoint2(bounds[0], bounds[3], 0);
//    gridPoints->Update();
//
//
//    // Perform the interpolation
//    vtkSmartPointer<vtkProbeFilter> probeFilter =  vtkSmartPointer<vtkProbeFilter>::New();
//    probeFilter->SetSourceData(filt->GetOutput());
//#if VTK_MAJOR_VERSION <= 5
//    probeFilter->SetInput(gridPoints->GetOutput());
//#else
//    probeFilter->SetInputConnection(gridPoints->GetOutputPort());
//#endif
//    probeFilter->Update();
//
//    vtkSmartPointer<vtkXMLPolyDataWriter> gridWriter =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//    gridWriter->SetFileName ( file_name.c_str());
//#if VTK_MAJOR_VERSION <= 5
//    gridWriter->SetInput(probeFilter->GetOutput());
//#else
//    probeFilter->SetInputConnection(probeFilter->GetOutputPort());
//#endif
//
//    gridWriter->Write();

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

    for(auto v = domain->vertices_begin(); v!=domain->vertices_end();v++ )
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
