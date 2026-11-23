// SPDX-License-Identifier: GPL-2.0-or-later
#include "checkpoint.hpp"
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
#include <optional>

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
    bool vtu = false, resume = false;
};

void partition(const Options& options, int worker, int workers)
{
    std::vector<fs::path> inputs{options.mesh};
    inputs.insert(inputs.end(), options.params.begin(), options.params.end());
    pt::ptree recipe;
    if (worker == 0)
    {
        recipe.put("format", 1);
        recipe.add_child("inputs", checkpoint::files(inputs));
        recipe.put("target_ranks", options.ranks);
        recipe.put("standalone", options.standalone);
        recipe.put("max_ghost_distance", options.distance);
        recipe.put("vtu", options.vtu);
        if (options.vtu)
            recipe.put("vtu_directory", fs::current_path().string());
        pt::ptree kept;
        for (int r : options.selected)
            kept.push_back({"", pt::ptree(std::to_string(r))});
        recipe.add_child("selected", kept);
    }
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
    const fs::path state = dir / ("." + manifest.filename().string());
    std::optional<checkpoint::Lock> lock;
    std::string generation(36, ' ');
    if (worker == 0)
    {
        lock.emplace(dir);
        fs::create_directories(state);
        generation = checkpoint::prepare(state / "run.json", recipe, options.resume);
        if (generation.size() != 36)
            throw std::runtime_error("Invalid partition checkpoint generation");
        fs::remove(manifest);
        std::cout << "Partitioning " << n << " faces for " << ranks.size() << " output ranks with " << workers
                  << " MPI workers\n"
                  << std::flush;
        const std::vector<fs::path> overview{"partition.vtu"};
        if (options.vtu && !(options.resume && checkpoint::complete(state / "overview.json", generation, overview)))
        {
            write_vtu("partition.vtu.partial", input);
            fs::rename("partition.vtu.partial", overview.front());
            checkpoint::finish(state / "overview.json", generation, overview);
        }
    }
    MPI_Bcast(generation.data(), 36, MPI_CHAR, 0, MPI_COMM_WORLD);
    int reused = 0;
    for (size_t i = worker; i < ranks.size(); i += workers)
    {
        const int rank = ranks[i];
        auto base = (dir / basename(rank)).string();
        std::vector<fs::path> outputs{base + "_mesh.h5", base + "_param.h5"};
        if (options.vtu)
            outputs.emplace_back("rank." + std::to_string(rank) + ".vtu");
        auto marker = state / (std::to_string(rank) + ".done.json");
        if (options.resume && checkpoint::complete(marker, generation, outputs))
        {
            ++reused;
            continue;
        }
        fs::remove(marker);
        auto local = select_rank(input, rank, ranks, options.distance, standalone);
        write_mesh(outputs[0].string() + ".partial", input, local, sizes, standalone);
        parameters.write(outputs[1].string() + ".partial", local.rows);
        if (options.vtu)
            write_vtu(outputs[2].string() + ".partial", input, &local);
        for (const auto& path : outputs)
            fs::rename(path.string() + ".partial", path);
        checkpoint::finish(marker, generation, outputs);
    }
    int total_reused;
    MPI_Reduce(&reused, &total_reused, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (worker == 0)
    {
        if (recipe.get_child("inputs") != checkpoint::files(inputs))
            throw std::runtime_error("Inputs changed while partitioning. Start a fresh run.");
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
        checkpoint::publish(manifest, tree);
        std::cout << "Completed " << ranks.size() - total_reused << " partitions, reused " << total_reused << ". Wrote "
                  << manifest << '\n';
    }
}

int run(int argc, char** argv, int worker, int workers)
{
    Options options;
    po::options_description desc("CHM partition options");
    desc.add_options()("help,h", "Show options")("resume", po::bool_switch(&options.resume),
                                                 "Reuse finished HDF5 partitions from an interrupted run")(
        "mesh-file,m", po::value(&options.mesh)->required(), "Input mesh")(
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
        if (options.standalone < 0 || options.standalone >= options.ranks)
            throw std::runtime_error("--standalone requires one valid target rank");
        if (!options.selected.empty() &&
            std::find(options.selected.begin(), options.selected.end(), options.standalone) == options.selected.end())
            throw std::runtime_error("--standalone rank is outside --valid-ranks");
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
    if (json && options.resume)
        throw std::runtime_error("Use --resume with the converted _mesh.h5 and _param.h5 inputs");
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
