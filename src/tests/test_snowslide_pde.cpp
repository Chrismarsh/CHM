// SPDX-License-Identifier: GPL-3.0-or-later
#include "snowslide_pde.hpp"
#include <boost/property_tree/json_parser.hpp>
#include <gtest/gtest.h>

namespace
{
constexpr double rad = std::numbers::pi / 180;
snowslide::Release paper(double slope, double swe, double vert, double vegetation = 0)
{
    return snowslide::release(snowslide::Profile::bernhardt_schulz, slope * rad, vegetation, 3178.4, -1.998, swe,
                              vert * std::cos(slope * rad), vert, 500);
}

TEST(SnowSlidePhysics, PaperTriggerDepthVegetationAndSteepSweFloor)
{
    EXPECT_EQ(paper(25, 1000, 10).mass, 0);
    EXPECT_GT(paper(25.001, 1000, 10).mass, 0);
    EXPECT_NEAR(paper(40, 800, 4).mass, 440, 1e-10);
    EXPECT_NEAR(paper(40, 800, 4, 3).mass, 200, 1e-10);
    EXPECT_NEAR(paper(74.99, 200, 1).mass, 100, 0.1);
    EXPECT_EQ(paper(75, 200, 1).mass, 150);
    EXPECT_EQ(paper(80, 200, 2).mass, 150);
    EXPECT_EQ(paper(80, 200, 0.5).mass, 150);
    EXPECT_EQ(paper(80, 35, 0.5).mass, 0);
    EXPECT_EQ(paper(80, 0, 0).capacity_left, 50);
}

TEST(SnowSlidePhysics, ChmNormalDepthAndDensity)
{
    const std::array slopes{30., 45., 60., 75.}, expected{3.079, 1.118, 0.445, 0.148};
    for (size_t i = 0; i < slopes.size(); ++i)
        EXPECT_NEAR(snowslide::chm_depth(slopes[i] * rad, 0, 3178.4, -1.998), expected[i], 0.001);
    const auto r = snowslide::release(snowslide::Profile::chm_vionnet, 60 * rad, 1.0, 3178.4, -1.998, 800, 2, 4, 500);
    EXPECT_EQ(r.mass, 400);
    EXPECT_EQ(r.normal_removed, 1);
    EXPECT_EQ(r.vertical_left, 2);
    EXPECT_EQ(r.capacity_left, 0);
    EXPECT_EQ(snowslide::chm_depth(5 * rad, 0, 3178.4, -1.998), 20);
}

TEST(SnowSlidePhysics, FluxReversesAndJacobianIncludesTheMovingSurface)
{
    const double a = 40, b = 20, rho = 300, T = 7, epsilon = 1e-5;
    for (double bed : {-2., 2.})
    {
        auto flux = [&](double x, double y) { return snowslide::edge_flux(x, y, bed + (x - y) / rho, T, rho); };
        const auto f = flux(a, b);
        EXPECT_DOUBLE_EQ(f.value, -snowslide::edge_flux(b, a, -bed + (b - a) / rho, T, rho).value);
        EXPECT_NEAR(f.donor_derivative, (flux(a + epsilon, b).value - flux(a - epsilon, b).value) / (2 * epsilon),
                    1e-8);
        EXPECT_NEAR(f.neighbour_derivative, (flux(a, b + epsilon).value - flux(a, b - epsilon).value) / (2 * epsilon),
                    1e-8);
    }
    EXPECT_EQ(snowslide::edge_flux(a, b, 0, T, rho).value, 0);
}

TEST(SnowSlidePhysics, TrackHoldsOnlyItsRemainingCapacity)
{
    EXPECT_EQ(snowslide::deposit(100, 1, 0, true), 0);
    EXPECT_EQ(snowslide::deposit(100, 1, 12, true), 12);
    EXPECT_NEAR(snowslide::deposit(100, 1, 0, false), 100 * (1 - std::exp(-1.)), 1e-12);
}

TEST(SnowSlidePhysics, WindErosionReducesMassAndBothDepthsBeforeRelease)
{
    double swe = 200, normal = 0.4, vertical = 2;
    snowslide::reserve_wind_erosion(75, swe, normal, vertical);
    EXPECT_DOUBLE_EQ(swe, 200);
    EXPECT_DOUBLE_EQ(normal, 0.4);
    EXPECT_DOUBLE_EQ(vertical, 2);
    snowslide::reserve_wind_erosion(-100, swe, normal, vertical);
    EXPECT_DOUBLE_EQ(swe, 100);
    EXPECT_DOUBLE_EQ(normal, 0.2);
    EXPECT_DOUBLE_EQ(vertical, 1);
    EXPECT_DOUBLE_EQ(paper(80, swe, vertical).mass, 50);
    snowslide::reserve_wind_erosion(-150, swe, normal, vertical);
    EXPECT_DOUBLE_EQ(swe, 0);
    EXPECT_DOUBLE_EQ(normal, 0);
    EXPECT_DOUBLE_EQ(vertical, 0);
}

// The JSON reader loads a serial mesh. Give its four faces the same ownership,
// halo metadata and row numbering used by CHM's partitioned mesh reader.
class FanMesh : public triangulation
{
  public:
    void distribute()
    {
        partition_mesh();
        for (auto f : _faces)
        {
            f->owner = determine_owner_of_global_index(static_cast<int>(f->cell_global_id), _num_faces_in_partition);
            f->is_ghost = f->owner != _comm_world.rank();
            if (!f->is_ghost)
                _global_to_local_faces_index_map[f->cell_global_id] = f->cell_local_id;
        }
        determine_local_boundary_faces();
        determine_process_ghost_faces_nearest_neighbors();
        determine_ghost_owners();
        setup_nearest_neighbor_communication();
    }
};

enum class Boundary
{
    closed,
    closed_with_wind,
    outlet,
    raised_fan,
    water
};

class SnowSlideFan : public ::testing::TestWithParam<Boundary>
{
};

TEST_P(SnowSlideFan, TwoSourcesFeedTwoFanFacesAndCloseTheWaterBudget)
{
    auto fan = std::make_shared<FanMesh>();
    if (fan->_comm_world.size() > 4)
        GTEST_SKIP() << "The fixture has four faces";
    std::istringstream json(R"({"mesh":{"version":"1.0.0","proj4":"+proj=utm +zone=11","is_geographic":0,
        "nvertex":6,"nelem":4,"vertex":[[0,0,0],[30,0,0],[15,26,0],[-15,26,120],[45,26,120],[15,-26,-3]],
        "elem":[[0,2,3],[0,1,2],[2,1,4],[1,0,5]],"neigh":[[-1,-1,1],[2,0,3],[-1,-1,1],[-1,-1,1]]},
        "parameters":{"snowslide_outlet":[0,0,0,0]}})");
    pt::ptree spec;
    pt::read_json(json, spec);
    const bool outlet = GetParam() == Boundary::outlet;
    const bool water = GetParam() == Boundary::water;
    const bool wind = GetParam() == Boundary::closed_with_wind;
    const bool trace = GetParam() == Boundary::closed || wind;
    if (GetParam() == Boundary::raised_fan)
        spec.get_child("mesh.vertex").back().second.back().second.put_value(0.3); // A 0.1 m adverse bed step.
    if (water)
    {
        pt::ptree cover;
        for (int id : {0, 0, 0, 1})
        {
            pt::ptree value;
            value.put_value(id);
            cover.push_back({"", value});
        }
        spec.add_child("parameters.landcover", cover);
    }
    fan->from_json(spec);
    fan->face(3)->parameter("snowslide_outlet"_s) = outlet ? 1 : 0;
    mesh domain = fan;
    config_file config;
    if (wind)
        config.put("reserve_drift_erosion", true);
    snowslide_pde slide(config);
    slide.global_param = std::make_shared<global>();
    slide.global_param->parameters.put("landcover.1.is_water", true);
    std::set<std::string> variables{"swe", "snowdepthavg", "snowdepthavg_vert"};
    if (wind)
        variables.insert("drift_mass");
    if (trace)
    {
        slide.set_optional_found("crocus_pending_snow");
        variables.insert("crocus_pending_snow");
    }
    for (const auto& v : *slide.provides())
        variables.insert(v.name);
    fan->init_timeseries(variables);
    fan->distribute();
    slide.init(domain);
    for (size_t i = 0; i < fan->size_local_faces(); ++i)
    {
        auto f = fan->face(i);
        const bool source = f->cell_global_id == 0 || f->cell_global_id == 2;
        (*f)["swe"_s] = source ? 200 : 0;
        (*f)["snowdepthavg"_s] = source ? 0.4 : 0;
        (*f)["snowdepthavg_vert"_s] = (*f)["snowdepthavg"_s] / std::cos(f->slope());
        if (wind)
        {
            constexpr std::array drift{-100., 0., 100., -9999.};
            (*f)["drift_mass"_s] = drift[f->cell_global_id];
        }
        if (trace)
        {
            (*f)["crocus_pending_snow"_s] = f->cell_global_id == 1 ? 0.00005 : 0;
            (*f)["swe"_s] += (*f)["crocus_pending_snow"_s];
        }
    }
    if (water)
        for (size_t i = 0; i < fan->size_local_faces(); ++i)
            if (auto f = fan->face(i); f->cell_global_id == 3)
                (*f)["swe"_s] = (*f)["snowdepthavg"_s] = (*f)["snowdepthavg_vert"_s] = -9999.;
    slide.global_param->first_time_step = false;
    slide.run(domain);
    // A closed and an exporting event must give the same fan split on all partitions.
    std::vector<double> adaptive_deposit(fan->size_local_faces());
    double net = 0, exported = 0;
    for (size_t i = 0; i < fan->size_local_faces(); ++i)
    {
        auto f = fan->face(i);
        const bool source = f->cell_global_id == 0 || f->cell_global_id == 2;
        const double source_release = wind && f->cell_global_id == 0 ? 50. : 150.;
        EXPECT_DOUBLE_EQ((*f)["avalanche_release_swe"_s], source ? source_release : 0.);
        if (water && f->cell_global_id == 3)
            EXPECT_EQ((*f)["avalanche_deposit_swe"_s], 0.);
        else if (!source)
            EXPECT_GT((*f)["avalanche_deposit_swe"_s], 1.0);
        EXPECT_GE((*f)["avalanche_deposit_swe"_s], 0);
        if (trace && f->cell_global_id == 1)
            EXPECT_DOUBLE_EQ((*f)["swe"_s], 0.00005);
        net += (*f)["delta_avalanche_mass"_s];
        exported += (*f)["avalanche_export_volume"_s];
        RecordProperty("face_" + std::to_string(f->cell_global_id) + "_deposit_mm",
                       std::to_string((*f)["avalanche_deposit_swe"_s]));
        adaptive_deposit[i] = (*f)["avalanche_deposit_swe"_s];
    }
    net = boost::mpi::all_reduce(fan->_comm_world, net, std::plus<double>());
    exported = boost::mpi::all_reduce(fan->_comm_world, exported, std::plus<double>());
    EXPECT_NEAR(net + exported, 0, 1e-7);
    if (outlet)
        EXPECT_GT(exported, 0);
    else
        EXPECT_EQ(exported, 0);
    // Persist the coupling state before another module consumes it.
    const auto path =
        boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("chm-slide-%%%%-%%%%-%%%%.nc");
    constexpr std::array keys{"delta_avalanche_mass", "delta_avalanche_snowdepth", "delta_avalanche_mass_sum",
                              "delta_avalanche_snowdepth_sum"};
    std::vector<std::array<double, 4>> saved(fan->size_local_faces());
    {
        netcdf file;
        file.create(path.string());
        slide.checkpoint(domain, file);
    }
    for (size_t i = 0; i < saved.size(); ++i)
        for (size_t k = 0; k < keys.size(); ++k)
        {
            auto& value = (*fan->face(i))[keys[k]];
            saved[i][k] = value;
            value = 0;
        }
    {
        netcdf file;
        file.open(path.string());
        slide.load_checkpoint(domain, file);
    }
    boost::filesystem::remove(path);
    for (size_t i = 0; i < saved.size(); ++i)
        for (size_t k = 0; k < keys.size(); ++k)
            EXPECT_DOUBLE_EQ((*fan->face(i))[keys[k]], saved[i][k]);

    // Repeat the same snapshot at a much finer fixed event step.
    config_file fixed = config;
    fixed.put("adaptive", false);
    fixed.put("event_step", 0.005);
    snowslide_pde reference(fixed);
    if (trace)
        reference.set_optional_found("crocus_pending_snow");
    reference.global_param = slide.global_param;
    reference.init(domain);
    reference.run(domain);
    for (size_t i = 0; i < saved.size(); ++i)
    {
        auto f = fan->face(i);
        EXPECT_NEAR((*f)["avalanche_deposit_swe"_s], adaptive_deposit[i], 2.0);
        (*f)["swe"_s] = 0;
    }
    slide.run(domain);
    for (size_t i = 0; i < fan->size_local_faces(); ++i)
        EXPECT_EQ((*fan->face(i))["delta_avalanche_mass"_s], 0);
    if (trace)
    {
        // A positive SWE with zero depths must be accounted for by the column's pending snow.
        for (size_t i = 0; i < fan->size_local_faces(); ++i)
            if (auto f = fan->face(i); f->cell_global_id == 1)
                (*f)["swe"_s] = 1.;
        EXPECT_THROW(slide.run(domain), module_error);
    }
}
INSTANTIATE_TEST_SUITE_P(Boundary, SnowSlideFan,
                         ::testing::Values(Boundary::closed, Boundary::closed_with_wind, Boundary::outlet,
                                           Boundary::raised_fan, Boundary::water));
} // namespace
