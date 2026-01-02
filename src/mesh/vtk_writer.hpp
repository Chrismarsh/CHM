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

#pragma once

#include <map>
#include <string>
#include <vector>

#ifdef USE_SPARSEHASH
#include <sparsehash/dense_hash_map>
#endif

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkTriangle.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>

class triangulation;

class vtk_writer
{
public:
    vtk_writer(triangulation* mesh, bool include_vertex_global_id = false, bool include_elevation_only = false);

    void init_grid(const std::vector<std::string>& output_variables);
    void update_data(const std::vector<std::string>& output_variables);
    void write_vtu(const std::string& file_name);

private:
    triangulation* _mesh;
    bool _include_vertex_global_id;
    bool _include_elevation_only;

    vtkSmartPointer<vtkUnstructuredGrid> _vtk_unstructuredGrid;
    vtkSmartPointer<vtkUnsignedLongArray> _vtu_global_id;

#ifdef USE_SPARSEHASH
    google::dense_hash_map< std::string, vtkSmartPointer<vtkFloatArray>  > data;
    google::dense_hash_map< std::string, vtkSmartPointer<vtkFloatArray>  > vectors;
    google::dense_hash_map< std::string, vtkSmartPointer<vtkFloatArray>  > vertex_data;
#else
    std::map<std::string, vtkSmartPointer<vtkFloatArray> > data;
    std::map<std::string, vtkSmartPointer<vtkFloatArray> > vectors;
    std::map<std::string, vtkSmartPointer<vtkFloatArray> > vertex_data;
#endif
};
