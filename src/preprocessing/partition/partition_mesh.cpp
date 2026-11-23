// SPDX-License-Identifier: GPL-2.0-or-later
#include "partition_mesh.hpp"
#include "coordinates.hpp"
#include <CGAL/centroid.h>
#include <bit>
#include <cstdint>
#include <numeric>
#include <unordered_map>

namespace chm_partition
{
namespace
{
void unique_sorted(std::vector<int>& ids)
{
    std::sort(ids.begin(), ids.end());
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}

struct Face
{
    Triple elem, neighbors;
    math::gis::Point_3 center;
    int ghost = -1; // -1 denotes a probed face outside the output halo
};

class Geometry
{
  public:
    explicit Geometry(const MeshInput& mesh) : mesh(mesh) {}
    void load(std::vector<int> ids)
    {
        std::erase_if(ids, [&](int id) { return faces.contains(id); });
        unique_sorted(ids);
        if (ids.empty())
            return;
        auto elems = read_rows<Triple>(mesh.elems, ids);
        auto neighbors = read_rows<Triple>(mesh.neighbors, ids);
        std::vector<int> needed;
        for (const auto& tri : elems)
            for (int v : tri)
                if (!points.contains(v))
                    needed.push_back(v);
        unique_sorted(needed);
        if (!needed.empty() && (needed.front() < 0 || needed.back() >= mesh.vertex_count))
            throw std::runtime_error("Mesh contains an out-of-range vertex");
        auto coords = read_rows<XYZ>(mesh.vertices, needed);
        for (size_t i = 0; i < needed.size(); ++i)
            points.emplace(needed[i], coords[i]);

        auto point = [&](int v)
        {
            const auto& p = points.at(v);
            return math::gis::Point_3(p[0], p[1], p[2]);
        };
        for (size_t i = 0; i < ids.size(); ++i)
        {
            const auto& e = elems[i];
            for (int n : neighbors[i])
                if (n < -1 || n >= mesh.offsets.back() || n == ids[i])
                    throw std::runtime_error("Invalid neighbor of face " + std::to_string(ids[i]));
            faces.emplace(ids[i], Face{e, neighbors[i], CGAL::centroid(point(e[0]), point(e[1]), point(e[2]))});
        }
    }

    // Each bit is an independent boundary seed. A face reached from seed A can only expand
    // seed A's search if its centroid is inside A's radius, even when other searches overlap.
    void halo(const std::vector<int>& owned, double radius)
    {
        std::vector<int> seeds, nearest;
        for (int id : owned)
        {
            bool boundary = false;
            for (int n : faces.at(id).neighbors)
                if (n == -1 || !faces.contains(n))
                {
                    boundary = true;
                    if (n != -1)
                        nearest.push_back(n);
                }
            if (boundary)
                seeds.push_back(id);
        }
        load(nearest);
        for (int id : nearest)
            faces.at(id).ghost = 1;

        auto distance = mesh.geographic ? math::gis::distance_latlong : math::gis::distance_UTM;
        using Masks = std::unordered_map<int, uint64_t>;
        constexpr int bits = std::numeric_limits<uint64_t>::digits;
        for (size_t start = 0; start < seeds.size(); start += bits)
        {
            Masks seen, frontier;
            for (size_t j = 0; j < std::min<size_t>(bits, seeds.size() - start); ++j)
                frontier[seeds[start + j]] = uint64_t{1} << j;
            while (!frontier.empty())
            {
                std::vector<int> missing;
                for (auto [id, mask] : frontier)
                    if (!faces.contains(id))
                        missing.push_back(id);
                load(std::move(missing));
                Masks next;
                for (auto [id, mask] : frontier)
                {
                    mask &= ~seen[id];
                    seen[id] |= mask;
                    auto& face = faces.at(id);
                    uint64_t accepted = 0;
                    while (mask)
                    {
                        auto j = std::countr_zero(mask);
                        if (distance(faces.at(seeds[start + j]).center, face.center) <= radius)
                            accepted |= uint64_t{1} << j;
                        mask &= mask - 1;
                    }
                    if (!accepted)
                        continue;
                    if (face.ghost < 0)
                        face.ghost = 2;
                    for (int n : face.neighbors)
                        if (n != -1 && (accepted & ~seen[n]))
                            next[n] |= accepted & ~seen[n];
                }
                frontier.swap(next);
            }
        }
    }
    const MeshInput& mesh;
    std::unordered_map<int, Face> faces;
    std::unordered_map<int, XYZ> points;
};
} // namespace

MeshInput::MeshInput(const fs::path& path, int ranks)
    : file(path.string(), H5F_ACC_RDONLY), elems(file.openDataSet("/mesh/elem")),
      neighbors(file.openDataSet("/mesh/neighbor")), vertices(file.openDataSet("/mesh/vertex")),
      ids(file.openDataSet("/mesh/cell_global_id")), projection(string_attribute(file, "/mesh/proj4")),
      method(string_attribute(file, "/mesh/partition_method")), vertex_count(row_count(vertices))
{
    hbool_t flag{};
    file.openAttribute("/mesh/is_geographic").read(H5::PredType::NATIVE_HBOOL, &flag);
    geographic = flag;
    if (std::stoi(string_attribute(file, "/mesh/version")) < 2)
        throw std::runtime_error("Mesh version 2 or newer required");
    if (file.attrExists("/mesh/is_partition"))
    {
        file.openAttribute("/mesh/is_partition").read(H5::PredType::NATIVE_HBOOL, &flag);
        if (flag)
            throw std::runtime_error("Expected a global mesh; input is already partitioned");
    }
    auto ds = file.openDataSet("/mesh/local_sizes");
    if (row_count(ds) != ranks)
        throw std::runtime_error("Requested rank count differs from /mesh/local_sizes");
    sizes.resize(ranks);
    read_range(ds, 0, std::span(sizes));
    offsets.push_back(0);
    for (int size : sizes)
    {
        if (size <= 0 || size > std::numeric_limits<int>::max() - offsets.back())
            throw std::runtime_error("Invalid partition size in /mesh/local_sizes");
        offsets.push_back(offsets.back() + size);
    }
    if (row_count(elems) != offsets.back() || row_count(neighbors) != offsets.back() ||
        row_count(ids) != offsets.back())
        throw std::runtime_error("Mesh dataset lengths disagree with /mesh/local_sizes");
}

bool MeshInput::ordered_range(int begin, int end) const
{
    std::vector<int> block;
    while (begin < end)
    {
        block.resize(std::min(batch_size, end - begin));
        read_range(ids, begin, std::span(block));
        for (int id : block)
            if (id != begin++)
                return false;
    }
    return true;
}

int MeshInput::owner_at(int position) const
{
    return static_cast<int>(std::upper_bound(offsets.begin(), offsets.end(), position) - offsets.begin()) - 1;
}

RankMesh select_rank(const MeshInput& input, int rank, std::span<const int> ranks, double radius, bool standalone)
{
    std::vector<int> owned(input.sizes.at(rank));
    if (input.ordered)
        std::iota(owned.begin(), owned.end(), input.offsets[rank]);
    else
        read_range(input.ids, input.offsets[rank], std::span(owned));
    Geometry geo(input);
    for (int id : owned)
        if (id < 0 || id >= input.offsets.back())
            throw std::runtime_error("Invalid owned face ID");
    geo.load(owned);
    if (geo.faces.size() != owned.size())
        throw std::runtime_error("Duplicate owned face ID");
    for (int id : owned)
        geo.faces.at(id).ghost = 0;
    if (!standalone)
        geo.halo(owned, radius);

    // Older global meshes may store a permutation in cell_global_id. Resolve only this
    // partition's selected IDs, in one bounded scan, without allocating a global inverse.
    std::unordered_map<int, int> positions;
    if (!input.ordered)
    {
        std::vector<int> block;
        for (int off = 0; off < input.offsets.back();)
        {
            block.resize(std::min(batch_size, input.offsets.back() - off));
            read_range(input.ids, off, std::span(block));
            for (int id : block)
            {
                auto face = geo.faces.find(id);
                if (face != geo.faces.end() && face->second.ghost >= 0 && !positions.emplace(id, off).second)
                    throw std::runtime_error("Repeated face in /mesh/cell_global_id");
                ++off;
            }
        }
    }

    std::vector<int> proxy(input.sizes.size(), -1), base(ranks.size() + 1, 0);
    for (size_t i = 0; i < ranks.size(); ++i)
    {
        proxy.at(ranks[i]) = i;
        base[i + 1] = base[i] + input.sizes[ranks[i]];
    }
    struct Selected
    {
        int id, row, owner, local;
    };
    std::vector<Selected> selected;
    std::unordered_map<int, int> output_id;
    for (const auto& [row, f] : geo.faces)
    {
        if (f.ghost < 0)
            continue;
        int pos = input.ordered ? row : positions.at(row);
        int owner = input.owner_at(pos);
        int p = proxy[owner];
        if (p < 0)
            continue;
        int local = pos - input.offsets[owner];
        int id = ranks.size() == input.sizes.size() ? row : base[p] + local;
        selected.push_back({id, row, p, local});
        output_id.emplace(row, id);
    }
    std::sort(selected.begin(), selected.end(), [](const auto& a, const auto& b) { return a.id < b.id; });
    RankMesh result;
    for (const auto& s : selected)
    {
        const auto& f = geo.faces.at(s.row);
        result.rows.push_back(s.row);
        result.ids.push_back(s.id);
        result.owners.push_back(s.owner);
        result.ghosts.push_back(f.ghost);
        result.local_ids.push_back(f.ghost ? s.row : s.local);
        result.elems.push_back(f.elem);
        result.vertex_ids.insert(result.vertex_ids.end(), f.elem.begin(), f.elem.end());
        auto n = f.neighbors;
        for (int& id : n)
        {
            auto it = output_id.find(id);
            id = it == output_id.end() ? -1 : it->second;
        }
        result.neighbors.push_back(n);
    }
    unique_sorted(result.vertex_ids);
    std::unordered_map<int, int> vertex_index;
    for (int id : result.vertex_ids)
    {
        vertex_index.emplace(id, result.vertices.size());
        result.vertices.push_back(geo.points.at(id));
    }
    for (auto& e : result.elems)
        for (int& v : e)
            v = vertex_index.at(v);
    return result;
}
} // namespace chm_partition
