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

#include <boost/filesystem/path.hpp>
#include <boost/property_tree/ptree.hpp>
namespace pt = boost::property_tree;

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

#include "triangulation.hpp"

class triangulation;

class vtk_writer
{
public:
    vtk_writer(mesh m, bool include_vertex_global_id = false, bool include_elevation_only = false);

    static void init_pvd(pt::ptree& pvd);
    static void append_pvd_entry(pt::ptree& pvd,
                                 const boost::filesystem::path& output_folder_path,
                                 const std::string& base_name,
                                 int rank,
                                 long timestep);
    static void write_pvd(const pt::ptree& pvd,
                          const boost::filesystem::path& output_folder_path,
                          const std::string& base_name);

    void set_write_ghost_neighbors(bool write_ghost_neighbors);
    void init_grid(const std::vector<std::string>& output_variables);
    void update_data(const std::vector<std::string>& output_variables);
    void write_vtu(const std::string& file_name);

private:
    mesh _mesh;
    bool _include_vertex_global_id;
    bool _include_elevation_only;
    bool _write_ghost_neighbors;

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
