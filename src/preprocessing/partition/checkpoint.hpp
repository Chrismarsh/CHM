// SPDX-License-Identifier: GPL-2.0-or-later
#pragma once
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <cerrno>
#include <filesystem>
#include <fstream>
#include <sys/stat.h>
#include <vector>

namespace chm_partition
{
namespace checkpoint
{
namespace fs = std::filesystem;
namespace pt = boost::property_tree;

inline pt::ptree identity(const fs::path& path)
{
    pt::ptree result;
    result.put("path", fs::canonical(path).string());
    struct stat metadata{};
    if (::stat(path.c_str(), &metadata) != 0)
        throw fs::filesystem_error("stat", path, std::error_code(errno, std::generic_category()));
    result.put("size", metadata.st_size);
    result.put("inode", metadata.st_ino);
#ifdef __APPLE__
    auto changed = metadata.st_ctimespec;
#else
    auto changed = metadata.st_ctim;
#endif
    result.put("changed_seconds", changed.tv_sec);
    result.put("changed_nanoseconds", changed.tv_nsec);
    return result;
}

inline pt::ptree files(const std::vector<fs::path>& paths)
{
    pt::ptree list;
    for (const auto& path : paths)
        list.push_back({"", identity(path)});
    return list;
}

inline void publish(const fs::path& path, const pt::ptree& tree)
{
    fs::path staging = path.string() + ".partial";
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(staging);
    pt::write_json(out, tree);
    out.close();
    fs::rename(staging, path);
}

// Only worker zero owns the process lock. The kernel releases it if the job is killed.
// Keep its inode in place across runs so two jobs cannot lock different files at the same path.
class Lock
{
  public:
    explicit Lock(const fs::path& directory)
    {
        fs::create_directories(directory);
        const auto path = (directory / ".partition.lock").string();
        std::ofstream(path, std::ios::app).close();
        lock = boost::interprocess::file_lock(path.c_str());
        if (!lock.try_lock())
            throw std::runtime_error("Another partition job is writing " + directory.string());
    }

  private:
    boost::interprocess::file_lock lock;
};

inline std::string prepare(const fs::path& path, const pt::ptree& recipe, bool resume)
{
    pt::ptree checkpoint;
    if (resume)
    {
        if (!fs::exists(path))
            throw std::runtime_error("No partition checkpoint at " + path.string());
        pt::read_json(path.string(), checkpoint);
        if (checkpoint.get_child("recipe") != recipe)
            throw std::runtime_error("Partition inputs or options changed. Start a fresh run without --resume.");
        return checkpoint.get<std::string>("generation");
    }
    auto id = boost::uuids::to_string(boost::uuids::random_generator()());
    checkpoint.put("generation", id);
    checkpoint.add_child("recipe", recipe);
    publish(path, checkpoint);
    return id;
}

inline bool complete(const fs::path& marker, const std::string& generation, const std::vector<fs::path>& outputs)
{
    if (!fs::exists(marker))
        return false;
    try
    {
        pt::ptree tree;
        pt::read_json(marker.string(), tree);
        return tree.get<std::string>("generation") == generation && tree.get_child("outputs") == files(outputs);
    }
    catch (const pt::ptree_error&)
    {
        return false;
    }
    catch (const fs::filesystem_error&)
    {
        return false;
    }
}

inline void finish(const fs::path& marker, const std::string& generation, const std::vector<fs::path>& outputs)
{
    pt::ptree tree;
    tree.put("generation", generation);
    tree.add_child("outputs", files(outputs));
    publish(marker, tree);
}
} // namespace checkpoint
} // namespace chm_partition
