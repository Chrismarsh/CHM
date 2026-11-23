// SPDX-License-Identifier: GPL-2.0-or-later
#pragma once
#include "partition_mesh.hpp"
namespace chm_partition
{
void write_vtu(const fs::path& path, const MeshInput& mesh, const RankMesh* rank = nullptr);
}
