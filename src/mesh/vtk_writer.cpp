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

#include "vtk_writer.hpp"

#include <cmath>

#include "triangulation.hpp"

vtk_writer::vtk_writer(triangulation* mesh, bool include_vertex_global_id, bool include_elevation_only)
    : _mesh(mesh)
    , _include_vertex_global_id(include_vertex_global_id)
    , _include_elevation_only(include_elevation_only)
    , _vtk_unstructuredGrid(nullptr)
    , _vtu_global_id(nullptr)
{
#ifdef USE_SPARSEHASH
    data.set_empty_key("");
    vectors.set_empty_key("");
    vertex_data.set_empty_key("");
#endif
}

void vtk_writer::init_grid(const std::vector<std::string>& output_variables)
{
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
    if (_mesh->_write_ghost_neighbors_to_vtu)
    {
        triangles->Allocate(_mesh->size_local_faces() + _mesh->_ghost_faces.size());
    }
    else
    {
        triangles->Allocate(_mesh->size_local_faces());
    }

    vtkSmartPointer<vtkStringArray> proj4 = vtkSmartPointer<vtkStringArray>::New();
    proj4->SetNumberOfComponents(1);
    proj4->SetName("proj4");
    proj4->InsertNextValue(_mesh->_srs_wkt);

    double scale = _mesh->is_geographic() == true ? 100000. : 1.;

    std::map<int, int> global_to_local_vertex_id;
    std::vector<int> global_vertex_id;

    // npoints holds the total number of points
    int npoints = 0;
    for (size_t i = 0; i < _mesh->size_local_faces(); i++)
    {
        mesh_elem fit = _mesh->face(i);

        vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();

        // loop over vertices of a face
        for (int j = 0; j < 3; ++j)
        {
            auto vit = fit->vertex(j);
            int global_id = vit->get_id();
            // If point hasn't been seen yet, account for it
            if (global_to_local_vertex_id.find(global_id) == global_to_local_vertex_id.end())
            {
                global_to_local_vertex_id[global_id] = npoints;
                npoints++;
                points->InsertNextPoint(vit->point().x() * scale, vit->point().y() * scale, vit->point().z());
                if (_include_vertex_global_id)
                {
                    global_vertex_id.push_back(global_id);
                }
            }
            tri->GetPointIds()->SetId(j, global_to_local_vertex_id[global_id]);
        }

        triangles->InsertNextCell(tri);
    }

    if (_mesh->_write_ghost_neighbors_to_vtu)
    {
        /* Ghost neighbors */
        for (size_t i = 0; i < _mesh->_ghost_faces.size(); i++)
        {
            mesh_elem fit = _mesh->_ghost_faces[i];

            vtkSmartPointer<vtkTriangle> tri = vtkSmartPointer<vtkTriangle>::New();

            // loop over vertices of a face
            for (int j = 0; j < 3; ++j)
            {
                auto vit = fit->vertex(j);
                int global_id = vit->get_id();
                // If point hasn't been seen yet, account for it
                if (global_to_local_vertex_id.find(global_id) == global_to_local_vertex_id.end())
                {
                    global_to_local_vertex_id[global_id] = npoints;
                    npoints++;
                    points->InsertNextPoint(vit->point().x() * scale, vit->point().y() * scale, vit->point().z());
                    if (_include_vertex_global_id)
                    {
                        global_vertex_id.push_back(global_id);
                    }
                }
                tri->GetPointIds()->SetId(j, global_to_local_vertex_id[global_id]);
            }

            triangles->InsertNextCell(tri);
        }
    }

    _vtk_unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    _vtk_unstructuredGrid->SetPoints(points);
    _vtk_unstructuredGrid->SetCells(VTK_TRIANGLE, triangles);
    _vtk_unstructuredGrid->GetFieldData()->AddArray(proj4);

    // assume that all the faces have the same number of variables and the same types of variables
    // by this point this should be a fair assumption
    auto variables = output_variables.size() == 0 ? _mesh->face(0)->variables() : output_variables;
    for (auto& v : variables)
    {
        data[v] = vtkSmartPointer<vtkFloatArray>::New();
        data[v]->SetName(v.c_str());
    }

    _vtu_global_id = vtkSmartPointer<vtkUnsignedLongArray>::New();
    _vtu_global_id->SetName("global_id");

    if (_mesh->_write_parameters)
    {
        auto params = _mesh->face(0)->parameters();
        for (auto& v : params)
        {
            data["[param] " + v] = vtkSmartPointer<vtkFloatArray>::New();
            data["[param] " + v]->SetName(("[param] " + v).c_str());
        }

        auto ics = _mesh->face(0)->initial_conditions();
        for (auto& v : ics)
        {
            data["[ic] " + v] = vtkSmartPointer<vtkFloatArray>::New();
            data["[ic] " + v]->SetName(("[ic] " + v).c_str());
        }

        // handle elevation/aspect/slope
        data["Elevation"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Elevation"]->SetName("Elevation");

        data["Slope"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Slope"]->SetName("Slope");

        data["Aspect"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Aspect"]->SetName("Aspect");

        data["Area"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Area"]->SetName("Area");

        data["is_ghost"] = vtkSmartPointer<vtkFloatArray>::New();
        data["is_ghost"]->SetName("is_ghost");

        data["ghost_type"] = vtkSmartPointer<vtkFloatArray>::New();
        data["ghost_type"]->SetName("ghost_type");

#ifdef USE_MPI
        data["owner"] = vtkSmartPointer<vtkFloatArray>::New();
        data["owner"]->SetName("owner");
#endif
    }
    else if (_include_elevation_only)
    {
        data["Elevation"] = vtkSmartPointer<vtkFloatArray>::New();
        data["Elevation"]->SetName("Elevation");
    }

    auto vec = _mesh->face(0)->vectors();
    for (auto& v : vec)
    {
        vectors[v] = vtkSmartPointer<vtkFloatArray>::New();
        vectors[v]->SetName(v.c_str());
        vectors[v]->SetNumberOfComponents(3);
    }

    if (_include_vertex_global_id)
    {
        vertex_data["global_id"] = vtkSmartPointer<vtkFloatArray>::New();
        vertex_data["global_id"]->SetName("global_id");
        for (int i = 0; i < npoints; ++i)
        {
            vertex_data["global_id"]->InsertTuple1(i, global_vertex_id[i]);
        }
    }
}

void vtk_writer::update_data(const std::vector<std::string>& output_variables)
{
    // if we haven't inited yet, do so.
    if (!_vtk_unstructuredGrid || _mesh->_terrain_deformed)
    {
        init_grid(output_variables);
    }

    auto variables = output_variables.size() == 0 ? _mesh->face(0)->variables() : output_variables;
    auto params = _mesh->face(0)->parameters();
    auto ics = _mesh->face(0)->initial_conditions();
    auto vecs = _mesh->face(0)->vectors();

    for (size_t i = 0; i < _mesh->size_local_faces(); i++)
    {
        mesh_elem fit = _mesh->face(i);

        for (auto& v : variables)
        {
            double d = (*fit)[v];
            if (d == -9999.)
            {
                d = nan("");
            }

            data[v]->InsertTuple1(i, d);
        }

        // this is mandatory now
        _vtu_global_id->InsertTuple1(i, fit->cell_global_id);

        if (_mesh->_write_parameters)
        {
            for (auto& v : params)
            {
                double d = fit->parameter(v);
                if (d == -9999.)
                {
                    d = nan("");
                }
                data["[param] " + v]->InsertTuple1(i, d);
            }

            for (auto& v : ics)
            {
                double d = fit->get_initial_condition(v);
                if (d == -9999.)
                {
                    d = nan("");
                }
                data["[ic] " + v]->InsertTuple1(i, d);
            }

            data["Elevation"]->InsertTuple1(i, fit->get_z());
            data["Slope"]->InsertTuple1(i, fit->slope());
            data["Aspect"]->InsertTuple1(i, fit->aspect());
            data["Area"]->InsertTuple1(i, fit->get_area());
            data["is_ghost"]->InsertTuple1(i, fit->is_ghost);
            data["ghost_type"]->InsertTuple1(i, fit->ghost_type);

#ifdef USE_MPI
            data["owner"]->InsertTuple1(i, _mesh->_comm_world.rank());
#endif
        }
        else if (_include_elevation_only)
        {
            data["Elevation"]->InsertTuple1(i, fit->get_z());
        }

        for (auto& v : vecs)
        {
            Vector_3 d = fit->face_vector(v);

            vectors[v]->InsertTuple3(i, d.x(), d.y(), d.z());
        }
    }

    if (_mesh->_write_ghost_neighbors_to_vtu)
    {
        /* Ghost neighbors */
        for (size_t i = 0; i < _mesh->_ghost_faces.size(); i++)
        {
            mesh_elem fit = _mesh->_ghost_faces[i];

            size_t insert_offset = i + _mesh->size_local_faces();

            for (auto& v : variables)
            {
                double d = -9999;
                if (fit->ghost_type == triangulation::GHOST_TYPE::NEIGH)
                {
                    d = (*fit)[v];
                }

                if (d == -9999.)
                {
                    d = nan("");
                }

                data[v]->InsertTuple1(insert_offset, d);
            }

            _vtu_global_id->InsertTuple1(insert_offset, fit->cell_global_id);

            if (_mesh->_write_parameters)
            {
                for (auto& v : params)
                {
                    double d = fit->parameter(v);
                    if (d == -9999.)
                    {
                        d = nan("");
                    }
                    data["[param] " + v]->InsertTuple1(insert_offset, d);
                }

                for (auto& v : ics)
                {
                    double d = fit->get_initial_condition(v);
                    if (d == -9999.)
                    {
                        d = nan("");
                    }
                    data["[ic] " + v]->InsertTuple1(insert_offset, d);
                }

                data["Elevation"]->InsertTuple1(insert_offset, fit->get_z());
                data["Slope"]->InsertTuple1(insert_offset, fit->slope());
                data["Aspect"]->InsertTuple1(insert_offset, fit->aspect());
                data["Area"]->InsertTuple1(insert_offset, fit->get_area());
                data["is_ghost"]->InsertTuple1(insert_offset, fit->is_ghost);
                data["ghost_type"]->InsertTuple1(insert_offset, fit->ghost_type);

                data["owner"]->InsertTuple1(insert_offset, fit->owner);
            }
            else if (_include_elevation_only)
            {
                data["Elevation"]->InsertTuple1(insert_offset, fit->get_z());
            }

            for (auto& v : vecs)
            {
                Vector_3 d = fit->face_vector(v);

                vectors[v]->InsertTuple3(insert_offset, d.x(), d.y(), d.z());
            }
        }
    }

    _vtk_unstructuredGrid->GetCellData()->AddArray(_vtu_global_id);

    for (auto& m : vectors)
    {
        _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);
    }

    for (auto& m : data)
    {
        _vtk_unstructuredGrid->GetCellData()->AddArray(m.second);
    }

    for (auto& m : vertex_data)
    {
        _vtk_unstructuredGrid->GetPointData()->AddArray(m.second);
    }
}

void vtk_writer::write_vtu(const std::string& file_name)
{
    // this now needs to be called from outside these functions
    // update_data();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(file_name.c_str());
//    writer->SetCompressorType( vtkXMLUnstructuredGridWriter::CompressorType::ZLIB);
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(_vtk_unstructuredGrid);
#else
    writer->SetInputData(_vtk_unstructuredGrid);
#endif
    writer->Write();
}
