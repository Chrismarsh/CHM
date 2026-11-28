// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once

#include "crocus_bridge.hpp"
#include "module_base.hpp"
#include <unordered_map>

/// ECCC Crocus multilayer snowpack with a prescribed ground temperature.
class Crocus : public module_base
{
    REGISTER_MODULE_HPP(Crocus);

  public:
    explicit Crocus(config_file cfg);
    void init(mesh& domain) override;
    void run(mesh_elem& face) override;
    void checkpoint(mesh& domain, netcdf& file) override;
    void load_checkpoint(mesh& domain, netcdf& file) override;

  private:
    void publish(mesh_elem& face);
    size_t row(const mesh_elem& face) const;
    std::unordered_map<size_t, size_t> pruned_rows_;
    int layers_;
    crocus::Boundary ground_;
    // Rank-local slabs, indexed by owned face, then field, then layer (surface first).
    std::vector<float> state_, pending_;
    std::vector<double> cumulative_subl_;
    std::vector<crocus::Fluxes> fluxes_;
};
