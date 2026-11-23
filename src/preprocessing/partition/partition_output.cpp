// SPDX-License-Identifier: GPL-2.0-or-later
#include "partition_mesh.hpp"
namespace chm_partition
{
void write_mesh(const fs::path& path, const MeshInput& input, const RankMesh& mesh, const std::vector<int>& sizes,
                bool standalone)
{
    H5::H5File file(path.string(), H5F_ACC_TRUNC);
    {
        auto group = file.createGroup("/mesh");
        // Every rank file records the owned counts for all selected ranks. The per-face arrays below
        // also contain this rank's ghosts, so their lengths can exceed its entry in local_sizes.
        write_values(group, "local_sizes", sizes);
        write_values(group, "owner", mesh.owners);
        write_values(group, "cell_global_id", mesh.ids);
        write_values(group, "ghost_type", mesh.ghosts);
        write_values(group, "elem", mesh.elems);
        write_values(group, "neighbor", mesh.neighbors);
        write_values(group, "vertex", mesh.vertices);
        write_attribute(file, "/mesh/version", "3.0.0");
        write_attribute(file, "/mesh/partition_method", input.method);
        write_attribute(file, "/mesh/proj4", input.projection);
        write_flag(file, "/mesh/is_geographic", input.geographic);
        write_flag(file, "/mesh/is_partition", !standalone);
    }
    file.close();
}

Parameters::Parameters(const std::vector<std::string>& paths, int faces)
{
    for (const auto& path : paths)
    {
        H5::H5File file(path, H5F_ACC_RDONLY);
        auto group = file.openGroup("/parameters");
        const auto count = group.getNumObjs();
        for (hsize_t i = 0; i < count; ++i)
        {
            auto name = group.getObjnameByIdx(i);
            auto ds = group.openDataSet(name);
            if (row_count(ds) != faces)
                throw std::runtime_error("Wrong parameter length: " + path + "/" + name);
            // Later parameter files override earlier definitions of the same field. Retain the source
            // group so each dataset and its HDF5 chunk cache can close after its column is copied.
            fields.insert_or_assign(name, group);
        }
    }
}
void Parameters::write(const fs::path& path, std::span<const int> rows) const
{
    H5::H5File file(path.string(), H5F_ACC_TRUNC);
    {
        auto group = file.createGroup("/parameters");
        hsize_t count = rows.size();
        for (const auto& [name, source_group] : fields)
        {
            auto src = source_group.openDataSet(name);
            auto dst = group.createDataSet(name, H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, &count));
            for (size_t start = 0; start < rows.size();)
            {
                // rows maps output positions to source faces. Read those faces in this order, then
                // write a contiguous slice starting at the matching position in the rank's parameter file.
                auto block = rows.subspan(start, std::min<size_t>(batch_size, rows.size() - start));
                auto values = read_rows<double>(src, block);
                hsize_t offset = start, length = values.size();
                auto space = dst.getSpace();
                space.selectHyperslab(H5S_SELECT_SET, &length, &offset);
                dst.write(values.data(), H5::PredType::NATIVE_DOUBLE, H5::DataSpace(1, &length), space);
                start += block.size();
            }
        }
    }
    file.close();
}
} // namespace chm_partition
