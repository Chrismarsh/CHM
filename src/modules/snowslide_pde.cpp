// SPDX-License-Identifier: GPL-3.0-or-later
#include "snowslide_pde.hpp"

REGISTER_MODULE_CPP(snowslide_pde);

namespace
{
// Volume diagnostics consumed by Snobal/FSM, including the pending increment on restart.
constexpr std::array coupling = {"delta_avalanche_mass", "delta_avalanche_snowdepth", "delta_avalanche_mass_sum",
                                 "delta_avalanche_snowdepth_sum"};
constexpr auto mobile_key = "snowslide_pde_mobile"_s;
constexpr auto head_key = "snowslide_pde_head"_s;
constexpr auto gx_key = "snowslide_pde_gradient_x"_s;
constexpr auto gy_key = "snowslide_pde_gradient_y"_s;

double global_max(mesh& domain, double value)
{
#ifdef USE_MPI
    return boost::mpi::all_reduce(domain->_comm_world, value, boost::mpi::maximum<double>());
#else
    return value;
#endif
}

double global_sum(mesh& domain, double value)
{
#ifdef USE_MPI
    return boost::mpi::all_reduce(domain->_comm_world, value, std::plus<double>());
#else
    return value;
#endif
}

struct XY
{
    double x, y;
};

template <class P, class Q> XY displacement(P a, Q b, bool geographic, double latitude)
{
    double dx = b.x() - a.x();
    double dy = b.y() - a.y();
    if (geographic)
    {
        constexpr double rad = std::numbers::pi / 180;
        constexpr double metres = 6378137.0 * rad; // same radius as math::gis::distance_latlong
        dx = std::remainder(dx, 360.0) * metres * std::cos(latitude * rad);
        dy *= metres;
    }
    return {dx, dy};
}
} // namespace

snowslide_pde::snowslide_pde(config_file cfg) : module_base("snowslide_pde", parallel::domain, cfg)
{
    const auto name = cfg.get<std::string>("release_profile", "bernhardt_schulz_2010");
    if (name == "bernhardt_schulz_2010")
        profile = snowslide::Profile::bernhardt_schulz;
    else if (name == "chm_vionnet_2021")
        profile = snowslide::Profile::chm_vionnet;
    else
        CHM_THROW_EXCEPTION(module_error, "Unknown snowslide_pde release_profile: " + name);

    mobility = cfg.get("mobility_length_m", 150.0);
    mobile_density = cfg.get("mobile_density", 300.0);
    deposit_density = cfg.get("deposit_density", 500.0);
    multiplier = cfg.get("avalache_mult", 3178.4);
    exponent = cfg.get("avalache_pow", -1.998);
    max_step = cfg.get("event_step", 0.25);
    min_step = cfg.get("min_event_step", 1e-5);
    adaptive = cfg.get("adaptive", true);
    abs_tol = cfg.get("event_abs_tol", 0.1);
    rel_tol = cfg.get("event_rel_tol", 0.01);
    for (double v : {mobility, mobile_density, deposit_density, multiplier, max_step, min_step, abs_tol, rel_tol})
        if (!std::isfinite(v) || v <= 0)
            CHM_THROW_EXCEPTION(module_error, "snowslide_pde requires finite positive transport and step parameters");
    if (!std::isfinite(exponent) || min_step > max_step)
        CHM_THROW_EXCEPTION(module_error, "Invalid snowslide_pde exponent or event step bounds");

    conflicts("snow_slide");
    for (auto key : {"swe", "snowdepthavg", "snowdepthavg_vert"})
        depends(key, SpatialType::neighbor);
    for (auto key : coupling)
        provides(key);
    for (auto key :
         {"maxDepth", "avalanche_release_swe", "avalanche_deposit_swe", "avalanche_export_volume",
          "snowslide_pde_mobile", "snowslide_pde_head", "snowslide_pde_gradient_x", "snowslide_pde_gradient_y"})
        provides(key);
}

void snowslide_pde::init(mesh& domain)
{
    const bool geographic = domain->is_geographic();
    const size_t count = domain->size_local_faces();
    cells.resize(count);
    export_rate.resize(count);
    outflow.resize(count);
    residual_swe.resize(count);
    bool invalid = false;
    for (size_t i = 0; i < count; ++i)
    {
        auto f = domain->face(i);
        auto& c = cells[i];
        c.area = f->get_area();
        c.slope = f->slope();
        c.cosine = std::max(0.001, std::cos(c.slope));
        c.water = is_water(f);
        c.vegetation = 0.0;
        if (!c.water && profile == snowslide::Profile::bernhardt_schulz)
        {
            const auto landcover = f->has_parameter("landcover"_s)
                                       ? "landcover." + std::to_string(static_cast<int>(f->parameter("landcover"_s)))
                                       : std::string{};
            if (f->has_parameter("snow_holding_depth"_s))
                c.vegetation = f->parameter("snow_holding_depth"_s);
            else if (auto value = global_param->parameters.get_optional<double>(landcover + ".snow_holding_depth"))
                c.vegetation = *value;
            else if (f->has_parameter("CanopyHeight"_s))
                c.vegetation = f->parameter("CanopyHeight"_s);
            else if (!landcover.empty())
                c.vegetation = global_param->parameters.get(landcover + ".CanopyHeight", 0.0);
        }
        else if (!c.water && f->has_vegetation())
            c.vegetation = f->veg_attribute("CanopyHeight");
        c.outlet = f->has_parameter("snowslide_outlet"_s) && f->parameter("snowslide_outlet"_s) > 0;
        invalid |= !(c.area > 0 && std::isfinite(c.area) && std::isfinite(c.slope) && c.vegetation >= 0 &&
                     std::isfinite(c.vegetation)) ||
                   (geographic && !f->has_parameter("area"_s));
        for (int j = 0; j < 3; ++j)
        {
            auto n = f->neighbor(j);
            const auto [a, b] = f->edge_vertexes(j);
            const double lat = (a.y() + b.y()) * 0.5;
            const auto tangent = displacement(a, b, geographic, lat);
            const auto to_a = displacement(f->center(), a, geographic, lat);
            const auto to_b = displacement(f->center(), b, geographic, lat);
            const XY delta = n != nullptr ? displacement(f->center(), n->center(), geographic, lat)
                                          : XY{(to_a.x + to_b.x) * 0.5, (to_a.y + to_b.y) * 0.5};
            auto& e = c.edges[j];
            e.dx = delta.x;
            e.dy = delta.y;
            e.length = std::hypot(tangent.x, tangent.y);
            e.nx = tangent.y / e.length;
            e.ny = -tangent.x / e.length;
            if (e.nx * e.dx + e.ny * e.dy < 0)
            {
                e.nx = -e.nx;
                e.ny = -e.ny;
            }
            e.normal_distance = e.nx * e.dx + e.ny * e.dy;
            invalid |= !(e.normal_distance > 0 && std::isfinite(e.normal_distance));
        }
        for (auto key : coupling)
            (*f)[key] = 0.0;
    }
    if (global_max(domain, invalid) > 0)
        CHM_THROW_EXCEPTION(module_error, "Invalid snowslide_pde mesh geometry or holding depth");
    problem = std::make_unique<math::LinearAlgebra::NearestNeighborProblem>(domain, 1);
}

double snowslide_pde::residual(mesh& domain, const Field& initial, const std::vector<double>& trial, double step,
                               bool assemble)
{
    for (size_t i = 0; i < cells.size(); ++i)
    {
        auto f = domain->face(i);
        const auto& c = cells[i];
        // A negative halo value marks a no-flux water face.
        (*f)[mobile_key] = c.water ? -1.0 : trial[i];
        (*f)[head_key] = f->center().z() + c.source.vertical_left +
                         initial[i].deposited / (deposit_density * c.cosine) + trial[i] / mobile_density;
    }
    domain->ghost_neighbors_communicate_variable(mobile_key);
    domain->ghost_neighbors_communicate_variable(head_key);

    // Weighted least-squares gradient; a lone neighbour resolves just the along-centre component.
    for (size_t i = 0; i < cells.size(); ++i)
    {
        auto f = domain->face(i);
        double xx = 0, xy = 0, yy = 0, bx = 0, by = 0;
        for (int j = 0; j < 3; ++j)
        {
            auto n = f->neighbor(j);
            if (n == nullptr || (*n)[mobile_key] < 0)
                continue;
            const auto& e = cells[i].edges[j];
            const double w = 1.0 / (e.dx * e.dx + e.dy * e.dy);
            const double dh = (*n)[head_key] - (*f)[head_key];
            xx += w * e.dx * e.dx;
            xy += w * e.dx * e.dy;
            yy += w * e.dy * e.dy;
            bx += w * e.dx * dh;
            by += w * e.dy * dh;
        }
        const double det = xx * yy - xy * xy;
        const double trace = xx + yy;
        if (det > 1e-12 * trace * trace)
        {
            (*f)[gx_key] = (yy * bx - xy * by) / det;
            (*f)[gy_key] = (xx * by - xy * bx) / det;
        }
        else
        {
            (*f)[gx_key] = trace > 0 ? bx / trace : 0.0;
            (*f)[gy_key] = trace > 0 ? by / trace : 0.0;
        }
    }
    domain->ghost_neighbors_communicate_variable(gx_key);
    domain->ghost_neighbors_communicate_variable(gy_key);

    if (assemble)
        problem->zeroSystem();
    double worst = 0;
    for (size_t i = 0; i < cells.size(); ++i)
    {
        auto f = domain->face(i);
        const auto& c = cells[i];
        double net = 0, diagonal = c.area;
        export_rate[i] = 0;
        outflow[i] = false;
        for (int j = 0; j < 3; ++j)
        {
            const auto& e = c.edges[j];
            auto n = f->neighbor(j);
            if (c.water || (n != nullptr && (*n)[mobile_key] < 0))
                continue;
            if (n == nullptr)
            {
                if (c.outlet)
                {
                    const double speed =
                        mobility * e.length * std::max(0.0, -((*f)[gx_key] * e.nx + (*f)[gy_key] * e.ny));
                    export_rate[i] += speed * trial[i] / c.area;
                    diagonal += step * speed;
                    outflow[i] |= speed > 0;
                }
                continue;
            }
            const double gx = 0.5 * ((*f)[gx_key] + (*n)[gx_key]);
            const double gy = 0.5 * ((*f)[gy_key] + (*n)[gy_key]);
            const double drop = (*f)[head_key] - (*n)[head_key];
            const double skew = gx * (e.dx - e.normal_distance * e.nx) + gy * (e.dy - e.normal_distance * e.ny);
            // Bound the reconstructed normal slope by the resolved gradient magnitude.
            const double limit =
                e.normal_distance * std::max(std::hypot(gx, gy), std::abs(drop) / std::hypot(e.dx, e.dy));
            const double corrected = std::clamp(drop + skew, -limit, limit);
            const auto flux = snowslide::edge_flux(trial[i], (*n)[mobile_key], corrected,
                                                   mobility * e.length / e.normal_distance, mobile_density);
            net += flux.value;
            outflow[i] |= corrected > 0;
            diagonal += step * flux.donor_derivative;
            if (assemble)
                problem->matrixSumIntoGlobalValues(f->cell_global_id, n->cell_global_id,
                                                   step * flux.neighbour_derivative);
        }
        residual_swe[i] = (trial[i] - initial[i].mobile) + step * (net / c.area + export_rate[i]);
        worst = std::max(worst, std::isfinite(residual_swe[i]) ? std::abs(residual_swe[i])
                                                               : std::numeric_limits<double>::infinity());
        if (assemble)
        {
            problem->matrixReplaceGlobalValues(f->cell_global_id, f->cell_global_id, diagonal);
            problem->rhsSumIntoGlobalValue(f->cell_global_id, -c.area * residual_swe[i]);
        }
    }
    return global_max(domain, worst);
}

bool snowslide_pde::advance(mesh& domain, Field& state, double step)
{
    std::vector<double> trial(state.size()), candidate(state.size());
    for (size_t i = 0; i < state.size(); ++i)
        trial[i] = state[i].mobile;

    constexpr double newton_tol = 1e-9; // residual in kg/m2, independent of the event error tolerance
    for (int iteration = 0; iteration < 30; ++iteration)
    {
        const double norm = residual(domain, state, trial, step, true);
        if (norm <= newton_tol)
        {
            // Recover the accepted inventory from the common edge fluxes. This prevents Krylov
            // residuals from accumulating as a distributed mass source over many substeps.
            bool invalid = false;
            for (size_t i = 0; i < state.size(); ++i)
                invalid |= trial[i] - residual_swe[i] < -newton_tol;
            if (global_max(domain, invalid) > 0)
                return false;
            for (size_t i = 0; i < state.size(); ++i)
            {
                const double transported = std::max(0.0, trial[i] - residual_swe[i]);
                const double settled = snowslide::deposit(
                    transported, step, cells[i].source.capacity_left - state[i].deposited, outflow[i]);
                state[i].mobile = transported - settled;
                state[i].deposited += settled;
                state[i].exported += step * export_rate[i];
            }
            return true;
        }
        if (!std::isfinite(norm))
            return false;
        const auto result = problem->Solve(false);
        if (global_max(domain, !result.converged) > 0)
            return false;
        const auto update = problem->getSolutionView();
        bool accepted = false;
        for (double damping = 1.0; damping >= 1.0 / 4096.0; damping *= 0.5)
        {
            bool invalid = false;
            for (size_t i = 0; i < state.size(); ++i)
            {
                candidate[i] = trial[i] + damping * update[domain->face(i)->cell_local_id];
                invalid |= !std::isfinite(candidate[i]) || candidate[i] < -newton_tol;
                candidate[i] = std::max(0.0, candidate[i]);
            }
            if (global_max(domain, invalid) == 0 && residual(domain, state, candidate, step, false) < norm)
            {
                trial.swap(candidate);
                accepted = true;
                break;
            }
        }
        if (!accepted)
            return false;
    }
    return false;
}

void snowslide_pde::run(mesh& domain)
{
    Field state(cells.size());
    double released = 0;
    bool invalid = false;
    for (size_t i = 0; i < cells.size(); ++i)
    {
        auto f = domain->face(i);
        auto& c = cells[i];
        double w = (*f)["swe"_s], h = (*f)["snowdepthavg"_s], v = (*f)["snowdepthavg_vert"_s];
        // CHM may call redistribution before the snow column publishes its first state.
        if (c.water || w == 0 || (global_param->first_time_step && is_nan(w)))
            w = h = v = 0;
        invalid |=
            !(w >= 0 && std::isfinite(w) && std::isfinite(h) && std::isfinite(v) && (w == 0 || (h > 0 && v > 0)));
        c.source = snowslide::release(profile, c.slope, c.vegetation, multiplier, exponent, w, h, v, deposit_density);
        state[i].mobile = c.source.mass;
        released += c.area * c.source.mass;
    }
    if (global_max(domain, invalid) > 0)
        CHM_THROW_EXCEPTION(module_error, "snowslide_pde received inconsistent SWE or snow depth");
    released = global_sum(domain, released);
    double remaining = released, step = max_step;
    size_t attempts = 0;
    while (remaining > released * 1e-8)
    {
        if (++attempts > 100000)
            CHM_THROW_EXCEPTION(module_error, "snowslide_pde event did not settle");
        Field next = state;
        bool ok = advance(domain, next, step);
        double error = 0;
        if (ok && adaptive)
        {
            Field fine = state;
            ok = advance(domain, fine, step * 0.5) && advance(domain, fine, step * 0.5);
            if (ok)
            {
                for (size_t i = 0; i < state.size(); ++i)
                {
                    auto scaled = [&](double a, double b)
                    { return std::abs(a - b) / (abs_tol + rel_tol * std::max(std::abs(a), std::abs(b))); };
                    error = std::max(
                        {error, scaled(next[i].mobile, fine[i].mobile),
                         scaled(next[i].deposited - state[i].deposited, fine[i].deposited - state[i].deposited),
                         scaled(next[i].exported - state[i].exported, fine[i].exported - state[i].exported)});
                }
                error = global_max(domain, error);
                next.swap(fine);
            }
        }
        if (!ok || error > 1.0)
        {
            if (step <= min_step)
                CHM_THROW_EXCEPTION(module_error,
                                    "snowslide_pde accuracy or Newton convergence failed at min_event_step");
            step = std::max(min_step, 0.5 * step);
            continue;
        }
        state.swap(next);
        remaining = 0;
        for (size_t i = 0; i < state.size(); ++i)
            remaining += cells[i].area * state[i].mobile;
        remaining = global_sum(domain, remaining);
        if (adaptive)
            step = std::clamp(step * std::clamp(0.9 / std::sqrt(std::max(1e-12, error)), 0.5, 2.0), min_step, max_step);
    }

    // Settle the tolerance-sized tail in place and check the event before publishing any delta.
    double accounted = 0;
    for (size_t i = 0; i < state.size(); ++i)
    {
        state[i].deposited += state[i].mobile;
        accounted += cells[i].area * (state[i].deposited + state[i].exported);
    }
    const double imbalance = global_sum(domain, accounted) - released;
    if (!std::isfinite(imbalance) || std::abs(imbalance) > 1e-7 * std::max(1.0, released))
        CHM_THROW_EXCEPTION(module_error, "snowslide_pde mass balance failed (kg): " + std::to_string(imbalance));
    for (size_t i = 0; i < state.size(); ++i)
    {
        auto f = domain->face(i);
        const auto& c = cells[i];
        const double mass = c.area * (state[i].deposited - c.source.mass) / 1000.0;
        const double depth = c.area * (state[i].deposited / deposit_density - c.source.normal_removed);
        (*f)["maxDepth"_s] = c.source.max_depth;
        (*f)["avalanche_release_swe"_s] = c.source.mass;
        (*f)["avalanche_deposit_swe"_s] = state[i].deposited;
        (*f)["avalanche_export_volume"_s] = c.area * state[i].exported / 1000.0;
        (*f)["delta_avalanche_mass"_s] = mass;
        (*f)["delta_avalanche_snowdepth"_s] = depth;
        (*f)["delta_avalanche_mass_sum"_s] += mass;
        (*f)["delta_avalanche_snowdepth_sum"_s] += depth;
    }
}

void snowslide_pde::checkpoint(mesh& domain, netcdf& file)
{
    for (auto key : coupling)
    {
        const std::string name = std::string("snowslide_pde:") + key;
        file.create_variable1D(name, domain->size_local_faces());
        for (size_t i = 0; i < domain->size_local_faces(); ++i)
            file.put_var1D(name, i, (*domain->face(i))[key]);
    }
}

void snowslide_pde::load_checkpoint(mesh& domain, netcdf& file)
{
    for (auto key : coupling)
    {
        const std::string name = std::string("snowslide_pde:") + key;
        for (size_t i = 0; i < domain->size_local_faces(); ++i)
            (*domain->face(i))[key] = file.get_var1D(name, i);
    }
}
