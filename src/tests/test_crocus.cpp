// SPDX-License-Identifier: GPL-3.0-or-later
#include "crocus.hpp"
#include <array>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>
#include <sstream>

TEST(Crocus, RestoresLayerStateAndPreviousStepDiagnostics)
{
    std::istringstream input(R"({"mesh":{"version":"1.0.0","nvertex":3,"nelem":1,"is_geographic":0,
        "proj4":"+proj=utm +zone=11 +datum=WGS84 +units=m +no_defs","vertex":[[600000,5600000,1000],
        [600100,5600000,1000],[600000,5600100,1000]],"elem":[[0,1,2]],"neigh":[[-1,-1,-1]]}})");
    pt::ptree description;
    pt::read_json(input, description);
    mesh domain = std::make_shared<triangulation>();
    domain->from_json(description);
    config_file cfg;
    cfg.put("max_layers", 3);
    cfg.put("ground_temperature_k", 268.0);
    Crocus model(cfg);
    model.global_param = std::make_shared<global>();
    std::set<std::string> fields;
    for (const auto& v : *model.provides())
        fields.insert(v.name);
    domain->init_timeseries(fields);
    // A point-mode selection can retain a nonzero mesh ID in its only row.
    domain->face(0)->cell_local_id = 37;
    model.init(domain);
    auto face = domain->face(0);
    EXPECT_EQ((*face)["Tsoil[0]"], 268.);
    EXPECT_EQ((*face)["swe"], 0.);
    const auto path = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("crocus-%%%%-%%%%.nc");
    struct Cleanup
    {
        boost::filesystem::path path;
        ~Cleanup()
        {
            boost::system::error_code ec;
            boost::filesystem::remove(path, ec);
        }
    } cleanup{path};
    constexpr int n = 3;
    std::array<float, crocus::field_count * n> state{};
    const float fluxes[] = {.73f, 12.f, -2.f, 260.f, .001f, .000003f, .0007f, -9.f};
    for (int k = 0; k < n; ++k)
    {
        state[k] = k + 1;
        state[crocus::density * n + k] = 200;
        state[crocus::heat * n + k] = -333700.f * state[k];
        state[crocus::diameter * n + k] = .0003f;
        state[crocus::sphericity * n + k] = .9f;
        state[crocus::history * n + k] = 1;
        state[crocus::age * n + k] = 10 + k;
    }
    {
        netcdf checkpoint;
        checkpoint.create(path.string());
        model.checkpoint(domain, checkpoint);
        auto& file = checkpoint.get_ncfile();
        file.getVar("crocus_state").putVar(state.data());
        file.getVar("crocus_fluxes").putVar(fluxes);
        const float trace = .0001f;
        const double cumulative = 7.25;
        file.getVar("crocus_pending_snow").putVar(&trace);
        file.getVar("crocus_cumulative_subl").putVar(&cumulative);
    }
    {
        netcdf checkpoint;
        checkpoint.open(path.string());
        model.load_checkpoint(domain, checkpoint);
        EXPECT_NEAR((*face)["swe"], 6.0001, 1e-7);
        EXPECT_EQ((*face)["crocus_pending_snow"], .0001f);
        EXPECT_NEAR((*face)["snowdepthavg"], .03, 1e-8);
        EXPECT_EQ((*face)["Nsnow"], 3);
        EXPECT_EQ((*face)["sum_snowpack_subl"], 7.25);
        EXPECT_EQ((*face)["snow_albedo"], fluxes[0]);
        EXPECT_EQ((*face)["roff"], fluxes[4]);
        cfg.put("ground_temperature_k", 263.f);
        Crocus other(cfg);
        EXPECT_THROW(other.load_checkpoint(domain, checkpoint), module_error);
    }
    {
        netcdf checkpoint;
        checkpoint.create(path.string());
        model.checkpoint(domain, checkpoint);
        std::array<float, crocus::field_count * n> restored{};
        checkpoint.get_ncfile().getVar("crocus_state").getVar(restored.data());
        EXPECT_EQ(restored, state);
    }
}
