// SPDX-License-Identifier: GPL-2.0-or-later
#pragma once

#include <H5Cpp.h>
#include <algorithm>
#include <array>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>

namespace chm_partition
{
using Triple = std::array<int, 3>;
using XYZ = std::array<double, 3>;
inline constexpr int batch_size = 8192;

inline const H5::DataType& memory_type(int*) { return H5::PredType::NATIVE_INT; }
inline const H5::DataType& memory_type(double*) { return H5::PredType::NATIVE_DOUBLE; }
template <class T> H5::ArrayType memory_type(std::array<T, 3>*)
{
    const hsize_t three = 3;
    return H5::ArrayType(memory_type(static_cast<T*>(nullptr)), 1, &three);
}

inline int row_count(const H5::DataSet& ds)
{
    hsize_t count;
    auto space = ds.getSpace();
    if (space.getSimpleExtentNdims() != 1 || space.getSimpleExtentDims(&count) != 1 ||
        count > std::numeric_limits<int>::max())
        throw std::runtime_error("Expected a one-dimensional CHM dataset with 32-bit row IDs: " + ds.getObjName());
    return static_cast<int>(count);
}

template <class T> void read_range(const H5::DataSet& ds, int start, std::span<T> values)
{
    if (values.empty())
        return;
    hsize_t offset = start, length = values.size();
    auto disk = ds.getSpace();
    disk.selectHyperslab(H5S_SELECT_SET, &length, &offset);
    ds.read(values.data(), memory_type(static_cast<T*>(nullptr)), H5::DataSpace(1, &length), disk);
}

// A contiguous owned run is a hyperslab. Sparse halo/vertex rows are one ordered point selection per batch.
template <class T> std::vector<T> read_rows(const H5::DataSet& ds, std::span<const int> rows)
{
    std::vector<T> values(rows.size());
    for (size_t i = 0; i < rows.size();)
    {
        auto block = rows.subspan(i, std::min<size_t>(batch_size, rows.size() - i));
        if (std::adjacent_find(block.begin(), block.end(), [](int a, int b) { return b - a != 1; }) == block.end())
            read_range(ds, block.front(), std::span(values).subspan(i, block.size()));
        else
        {
            std::vector<hsize_t> positions(block.begin(), block.end());
            hsize_t length = positions.size();
            auto disk = ds.getSpace();
            disk.selectElements(H5S_SELECT_SET, length, positions.data());
            ds.read(values.data() + i, memory_type(static_cast<T*>(nullptr)), H5::DataSpace(1, &length), disk);
        }
        i += block.size();
    }
    return values;
}

template <class T> void write_values(H5::Group& group, const char* name, const std::vector<T>& values)
{
    hsize_t n = values.size();
    const auto& native = memory_type(static_cast<T*>(nullptr));
    const H5::DataType& wire =
        std::is_same_v<T, int> ? static_cast<const H5::DataType&>(H5::PredType::STD_I32BE) : native;
    auto ds = group.createDataSet(name, wire, H5::DataSpace(1, &n));
    if (n)
        ds.write(values.data(), native);
}

inline std::string string_attribute(const H5::H5File& file, const char* name)
{
    auto attr = file.openAttribute(name);
    std::string value;
    attr.read(attr.getStrType(), value);
    return value;
}
inline void write_attribute(H5::H5File& file, const char* name, const std::string& value)
{
    const hsize_t one = 1;
    H5::StrType type(H5::PredType::C_S1, 256);
    if (value.size() >= type.getSize())
        throw std::runtime_error("CHM string attribute is too long: " + std::string(name));
    file.createAttribute(name, type, H5::DataSpace(1, &one)).write(type, value);
}
inline void write_flag(H5::H5File& file, const char* name, bool value)
{
    const hsize_t one = 1;
    hbool_t flag = value;
    file.createAttribute(name, H5::PredType::NATIVE_HBOOL, H5::DataSpace(1, &one))
        .write(H5::PredType::NATIVE_HBOOL, &flag);
}
} // namespace chm_partition
