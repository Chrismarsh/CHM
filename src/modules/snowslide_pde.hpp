// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include "LinearAlgebra.hpp"
#include "module_base.hpp"
#include "snowslide_physics.hpp"

/** \ingroup modules snow
 * Conservative finite-volume transport of mobile avalanche snow. Configuration, units,
 * and release-profile references are in docs/snowslide_pde.rst.
 */
class snowslide_pde : public module_base
{
    REGISTER_MODULE_HPP(snowslide_pde);

  public:
    explicit snowslide_pde(config_file cfg);
    void init(mesh& domain) override;
    void run(mesh& domain) override;
    void checkpoint(mesh& domain, netcdf& file) override;
    void load_checkpoint(mesh& domain, netcdf& file) override;

  private:
    struct Edge
    {
        double dx{}, dy{}, nx{}, ny{}, length{}, normal_distance{};
    };
    struct Cell
    {
        std::array<Edge, 3> edges;
        double area{}, slope{}, cosine{}, vegetation{};
        bool outlet{}, water{};
        snowslide::Release source{};
    };
    struct State
    {
        double mobile{}, deposited{}, exported{}; // kg/m2
    };
    using Field = std::vector<State>;
    // Evaluates the full nonlinear residual, including the reconstructed gradient.
    // Optionally builds the Jacobian of the nearest-neighbour part.
    double residual(mesh& domain, const Field& initial, const std::vector<double>& trial, double step, bool assemble);
    bool advance(mesh& domain, Field& state, double step);
    std::vector<Cell> cells;
    std::vector<double> export_rate;
    std::vector<double> residual_swe;
    std::vector<char> outflow;
    std::unique_ptr<math::LinearAlgebra::NearestNeighborProblem> problem;
    snowslide::Profile profile;
    double mobility, mobile_density, deposit_density, multiplier, exponent;
    double max_step, min_step, abs_tol, rel_tol;
    bool adaptive, reserve_drift_erosion;
};
