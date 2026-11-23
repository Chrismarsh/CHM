// SPDX-License-Identifier: GPL-2.0-or-later
#include "convert_mesh.hpp"
#include "partition_mesh.hpp"
#include "partition_vtu.hpp"

#include <boost/exception/diagnostic_information.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <numeric>

namespace
{
using namespace chm_partition;
namespace po = boost::program_options;
namespace pt = boost::property_tree;
struct Options
{
    std::string mesh;
    std::vector<std::string> params;
    std::vector<int> selected;
    int ranks = 0, standalone = -1;
    double distance = 100;
    bool vtu = false;
};

void partition(const Options& options, int worker, int workers)
{
    MeshInput input(options.mesh, options.ranks);
    Parameters parameters(options.params, input.offsets.back());
    // Distribute the one-time permutation check even when there are more workers than output files.
    int64_t n = input.offsets.back();
    int ordered = input.ordered_range(n * worker / workers, n * (worker + 1) / workers);
    int all_ordered;
    MPI_Allreduce(&ordered, &all_ordered, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    input.ordered = all_ordered;

    const auto& ranks = options.selected;
    std::vector<int> sizes;
    for (int rank : ranks)
        sizes.push_back(input.sizes[rank]);
    fs::path mesh(options.mesh);
    auto stem = mesh.stem().string();
    auto prefix = mesh.parent_path() / (stem + ".np" + std::to_string(options.ranks));
    fs::path dir = prefix.string() + ".partition.meshes";
    const bool standalone = options.standalone >= 0;
    fs::path manifest =
        prefix.string() + (standalone ? ".standalone." + std::to_string(options.standalone) : "") + ".partition";
    auto basename = [&](int rank)
    { return stem + ".partition." + (standalone ? "standalone." : "") + std::to_string(rank); };
    if (worker == 0)
    {
        fs::create_directories(dir);
        std::cout << "Partitioning " << n << " faces for " << ranks.size() << " output ranks with " << workers
                  << " MPI workers\n"
                  << std::flush;
        if (options.vtu)
            write_vtu("partition.vtu", input);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (size_t i = worker; i < ranks.size(); i += workers)
    {
        const int rank = ranks[i];
        auto local = select_rank(input, rank, ranks, options.distance, standalone);
        auto base = (dir / basename(rank)).string();
        write_mesh(base + "_mesh.h5", input, local, sizes, standalone);
        parameters.write(base + "_param.h5", local.rows);
        if (options.vtu)
            write_vtu("rank." + std::to_string(rank) + ".vtu", input, &local);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (worker == 0)
    {
        pt::ptree tree, meshes, params, local_sizes;
        tree.put("ranks", ranks.size());
        tree.put("max_ghost_distance", options.distance);
        tree.put("num_global_faces", std::accumulate(sizes.begin(), sizes.end(), int64_t{}));
        for (int rank : ranks)
        {
            auto base = (dir / basename(rank)).string();
            meshes.push_back({"", pt::ptree(base + "_mesh.h5")});
            params.push_back({"", pt::ptree(base + "_param.h5")});
            local_sizes.push_back({"", pt::ptree(std::to_string(input.sizes[rank]))});
        }
        tree.add_child("meshes", meshes);
        tree.add_child("parameters", params);
        tree.add_child("local_sizes", local_sizes);
        pt::write_json(manifest.string(), tree);
        std::cout << "Wrote " << manifest << '\n';
    }
}

int run(int argc, char** argv, int worker, int workers)
{
    Options options;
    po::options_description desc("CHM partition options");
    desc.add_options()("help,h", "Show options")("mesh-file,m", po::value(&options.mesh)->required(), "Input mesh")(
        "param-file,p", po::value(&options.params)->composing()->required(), "Input parameter file (repeatable)")(
        "mpi-ranks", po::value(&options.ranks), "Number of target partitions in the input mesh")(
        "max-ghost-distance,g", po::value(&options.distance), "Maximum ghost distance in metres (default 100)")(
        "valid-ranks,v", po::value(&options.selected)->composing(), "Output a subset of target ranks (repeatable)")(
        "standalone,s", po::value(&options.standalone), "Output an owned-only mesh for this rank")(
        "write-vtu", po::bool_switch(&options.vtu), "Write partition.vtu and rank.N.vtu");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if (vm.count("help"))
    {
        if (!worker)
            std::cout << desc << '\n';
        return 0;
    }
    po::notify(vm);
    bool json = fs::path(options.mesh).extension() == ".mesh";
    if (!json && fs::path(options.mesh).extension() != ".h5")
        throw std::runtime_error("Expected a .mesh or .h5 mesh file");
    if ((!json || vm.count("mpi-ranks")) && options.ranks < 2)
        throw std::runtime_error("--mpi-ranks must be at least 2");
    if (!std::isfinite(options.distance) || options.distance < 0)
        throw std::runtime_error("Invalid ghost distance");
    if (vm.count("standalone"))
    {
        if (!options.selected.empty() || options.standalone < 0 || options.standalone >= options.ranks)
            throw std::runtime_error(
                "--standalone requires one valid target rank and cannot be combined with --valid-ranks");
        options.selected = {options.standalone};
    }
    if (options.selected.empty())
    {
        options.selected.resize(options.ranks);
        std::iota(options.selected.begin(), options.selected.end(), 0);
    }
    else
    {
        auto& ranks = options.selected;
        std::sort(ranks.begin(), ranks.end());
        if (ranks.front() < 0 || ranks.back() >= options.ranks ||
            std::adjacent_find(ranks.begin(), ranks.end()) != ranks.end())
            throw std::runtime_error("Invalid or repeated target rank");
    }
    if (json)
    {
        if (worker == 0)
            convert_json_mesh(options.mesh, options.params);
        MPI_Barrier(MPI_COMM_WORLD);
        if (!options.ranks)
            return 0;
        auto stem = fs::path(options.mesh).stem().string();
        options.mesh = stem + "_mesh.h5";
        options.params = {stem + "_param.h5"};
    }
    partition(options, worker, workers);
    return 0;
}
} // namespace
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int worker, workers;
    MPI_Comm_rank(MPI_COMM_WORLD, &worker);
    MPI_Comm_size(MPI_COMM_WORLD, &workers);
    H5::Exception::dontPrint();
    try
    {
        run(argc, argv, worker, workers);
    }
    catch (const H5::Exception& e)
    {
        std::cerr << "partition worker " << worker << ": " << e.getDetailMsg() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    catch (...)
    {
        std::cerr << "partition worker " << worker << ": " << boost::current_exception_diagnostic_information() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
    MPI_Finalize();
    return 0;
}
