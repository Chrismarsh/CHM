// SPDX-License-Identifier: GPL-2.0-or-later
#pragma once
#include "hdf_rows.hpp"
#include <filesystem>
#include <map>

namespace chm_partition
{
namespace fs = std::filesystem;
struct MeshInput
{
    explicit MeshInput(const fs::path& path, int ranks);
    H5::H5File file;
    H5::DataSet elems, neighbors, vertices, ids;
    std::vector<int> sizes, offsets;
    std::string projection, method;
    bool geographic, ordered = true;
    int vertex_count;
    bool ordered_range(int begin, int end) const;
    int owner_at(int position) const;
};

struct RankMesh
{
    std::vector<int> rows, ids, owners, ghosts, vertex_ids, local_ids;
    std::vector<Triple> elems, neighbors;
    std::vector<XYZ> vertices;
};

RankMesh select_rank(const MeshInput& input, int rank, std::span<const int> ranks, double distance, bool standalone);
void write_mesh(const fs::path& path, const MeshInput& input, const RankMesh& mesh, const std::vector<int>& sizes,
                bool standalone);
class Parameters
{
  public:
    Parameters(const std::vector<std::string>& paths, int faces);
    void write(const fs::path& path, std::span<const int> rows) const;

  private:
    std::map<std::string, H5::DataSet> fields;
};
} // namespace chm_partition
