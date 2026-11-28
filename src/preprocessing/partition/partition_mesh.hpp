// SPDX-License-Identifier: GPL-2.0-or-later
#pragma once
#include "hdf_rows.hpp"
#include <filesystem>
#include <map>

namespace chm_partition
{
namespace fs = std::filesystem;
// Each worker keeps the input datasets open across output ranks. The ownership table has one size per target rank.
// Rank r owns positions [offsets[r], offsets[r + 1]) in /mesh/cell_global_id. The value at each position
// is the source row in elem, neighbor and every parameter column. For example, sizes {2, 1} and IDs
// {2, 0, 1} assign source rows 2 and 0 to rank 0. When ordered is true, position and row coincide.
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

// Face arrays share one output ordering, ascending by ids. rows[i] selects the original parameter values
// for output face i. ids[i] is the ID CHM uses to resolve its neighbors, including across output files.
// owners uses the selected ranks' new numbering. ghosts uses CHM's codes 0 (owned), 1 (direct), 2 (distance).
// For VTU diagnostics, local_ids holds the owner's local index for owned faces and the source row for ghosts.
struct RankMesh
{
    std::vector<int> rows, ids, owners, ghosts, vertex_ids, local_ids;
    // elems indexes the compact vertices array. vertex_ids maps those indices back to source vertices.
    // neighbors holds output face IDs, with -1 for an exterior edge or a neighbor absent from this file.
    std::vector<Triple> elems, neighbors;
    std::vector<XYZ> vertices;
};

// rank is the original METIS rank. ranks must contain it and lists the selected original ranks in ascending order.
RankMesh select_rank(const MeshInput& input, int rank, std::span<const int> ranks, double distance, bool standalone);
void write_mesh(const fs::path& path, const MeshInput& input, const RankMesh& mesh, const std::vector<int>& sizes,
                bool standalone);
class Parameters
{
  public:
    Parameters(const std::vector<std::string>& paths, int faces);
    void write(const fs::path& path, std::span<const int> rows) const;

  private:
    std::map<std::string, H5::Group> fields;
};
} // namespace chm_partition
