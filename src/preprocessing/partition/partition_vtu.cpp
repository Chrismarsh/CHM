// SPDX-License-Identifier: GPL-2.0-or-later
#include "partition_vtu.hpp"
#include <fstream>
#include <iomanip>
#include <numeric>

namespace chm_partition
{
// A null rank requests the global overview. Stream one VTK piece with shared vertex indices across
// batches. For a rank diagnostic, RankMesh already holds the compact connectivity and coordinates.
void write_vtu(const fs::path& path, const MeshInput& mesh, const RankMesh* rank)
{
    std::ofstream out;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out.open(path);
    out << std::setprecision(17);
    int nv = rank ? static_cast<int>(rank->vertices.size()) : mesh.vertex_count;
    int nf = rank ? static_cast<int>(rank->rows.size()) : mesh.offsets.back();
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"><UnstructuredGrid>\n"
        << "<FieldData><Array type=\"String\" Name=\"proj4\" NumberOfTuples=\"1\" format=\"ascii\">";
    for (unsigned char c : mesh.projection)
        out << static_cast<unsigned>(c) << ' ';
    out << "0</Array></FieldData>\n<Piece NumberOfPoints=\"" << nv << "\" NumberOfCells=\"" << nf << "\">\n";
    auto array = [&](const char* name, const char* type, int components, auto write)
    {
        out << "<DataArray Name=\"" << name << "\" type=\"" << type << "\" NumberOfComponents=\"" << components
            << "\" format=\"ascii\">\n";
        write();
        out << "\n</DataArray>\n";
    };
    auto integers = [&](const char* name, int size, auto value)
    {
        array(name, "Int64", 1,
              [&]
              {
                  for (int i = 0; i < size; ++i)
                      out << value(i) << ' ';
              });
    };
    out << "<Points>";
    array("Points", "Float64", 3,
          [&]
          {
              // Match CHM's vtk_writer display scale for lon/lat. HDF5 retains the source coordinates.
              double scale = mesh.geographic ? 100000. : 1.;
              for (int i = 0; i < nv;)
              {
                  int count = std::min(batch_size, nv - i);
                  std::vector<XYZ> block;
                  if (!rank)
                  {
                      block.resize(count);
                      read_range(mesh.vertices, i, std::span(block));
                  }
                  for (int j = 0; j < count; ++j)
                  {
                      const auto& p = rank ? rank->vertices[i + j] : block[j];
                      out << p[0] * scale << ' ' << p[1] * scale << ' ' << p[2] << '\n';
                  }
                  i += count;
              }
          });
    out << "</Points><PointData>";
    integers("global_id", nv, [&](int i) { return rank ? rank->vertex_ids[i] : i; });
    out << "</PointData><Cells>";
    // All overview cell arrays follow partition order. Use cell_global_id to fetch the corresponding
    // connectivity and elevations even when the source geometry rows are permuted.
    auto cells = [&](auto consume)
    {
        for (int i = 0; i < nf;)
        {
            int count = std::min(batch_size, nf - i);
            if (rank)
                consume(i, std::span(rank->elems).subspan(i, count));
            else
            {
                std::vector<int> rows(count);
                read_range(mesh.ids, i, std::span(rows));
                auto elems = read_rows<Triple>(mesh.elems, rows);
                consume(i, std::span<const Triple>(elems));
            }
            i += count;
        }
    };
    array("connectivity", "Int32", 1,
          [&]
          {
              cells(
                  [&](int, auto tris)
                  {
                      for (auto e : tris)
                          out << e[0] << ' ' << e[1] << ' ' << e[2] << '\n';
                  });
          });
    integers("offsets", nf, [](int i) { return 3LL * (i + 1LL); });
    array("types", "UInt8", 1,
          [&]
          {
              for (int i = 0; i < nf; ++i)
                  out << "5 ";
          });
    out << "</Cells><CellData>";
    integers("owner", nf, [&](int i) { return rank ? rank->owners[i] : mesh.owner_at(i); });
    array("global_id", "Int32", 1,
          [&]
          {
              if (rank)
              {
                  for (int id : rank->ids)
                      out << id << ' ';
              }
              else
                  for (int i = 0; i < nf;)
                  {
                      std::vector<int> rows(std::min(batch_size, nf - i));
                      read_range(mesh.ids, i, std::span(rows));
                      for (int row : rows)
                          out << row << ' ';
                      i += static_cast<int>(rows.size());
                  }
          });
    array("Elevation", "Float64", 1,
          [&]
          {
              cells(
                  [&](int, auto tris)
                  {
                      std::vector<int> ids;
                      std::vector<XYZ> points;
                      if (!rank)
                      {
                          for (auto e : tris)
                              ids.insert(ids.end(), e.begin(), e.end());
                          std::sort(ids.begin(), ids.end());
                          ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
                          points = read_rows<XYZ>(mesh.vertices, ids);
                      }
                      auto z = [&](int v)
                      {
                          return rank ? rank->vertices[v][2]
                                      : points[std::lower_bound(ids.begin(), ids.end(), v) - ids.begin()][2];
                      };
                      for (auto e : tris)
                          out << (z(e[0]) + z(e[1]) + z(e[2])) / 3. << ' ';
                  });
          });
    if (rank)
    {
        integers("local_id", nf, [&](int i) { return rank->local_ids[i]; });
        integers("ghost_type", nf, [&](int i) { return rank->ghosts[i]; });
        integers("is_ghost", nf, [&](int i) { return int(rank->ghosts[i] != 0); });
    }
    out << "</CellData></Piece></UnstructuredGrid></VTKFile>\n";
    out.close();
}
} // namespace chm_partition
