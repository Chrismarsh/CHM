// SPDX-License-Identifier: GPL-2.0-or-later
#include "convert_mesh.hpp"
#include "readjson.hpp"
#include "triangulation.hpp"
#include <filesystem>

namespace
{
class JsonMesh : public triangulation
{
  public:
    void convert(const std::string& path, const std::vector<std::string>& parameters)
    {
        auto mesh = read_json(path);
        _is_geographic = mesh.get<int>("mesh.is_geographic") == 1;
        math::gis::distance = _is_geographic ? math::gis::distance_latlong : math::gis::distance_UTM;
        math::gis::point_from_bearing =
            _is_geographic ? math::gis::point_from_bearing_latlong : math::gis::point_from_bearing_UTM;
        for (const auto& parameter : parameters)
        {
            auto values = read_json(parameter);
            for (const auto& [key, column] : values)
                mesh.put_child("parameters." + key, column);
        }
        if (_is_geographic && !mesh.get_child_optional("parameters.area"))
            throw std::runtime_error("Geographic meshes require the area parameter");
        if (!mesh.get_child_optional("mesh.cell_global_id"))
            SPDLOG_WARN("No face permutation was found in the mesh");
        from_json(mesh);
        _num_faces = _num_global_faces = _faces.size();
        _local_faces = _faces;
        to_hdf5(std::filesystem::path(path).stem().string());
    }
};
} // namespace
void convert_json_mesh(const std::string& path, const std::vector<std::string>& parameters)
{
    JsonMesh mesh;
    mesh.convert(path, parameters);
}
