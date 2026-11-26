// SPDX-License-Identifier: GPL-3.0-or-later
#include "crocus.hpp"

#include <array>
#include <cmath>
#include <meteoio/MeteoIO.h>
#include <mutex>
#include <numbers>

REGISTER_MODULE_CPP(Crocus);

namespace
{
constexpr std::array snow_temp = {"Tsnow[0]"_s, "Tsnow[1]"_s, "Tsnow[2]"_s};
constexpr std::array snow_liquid = {"Sliq[0]"_s, "Sliq[1]"_s, "Sliq[2]"_s};
constexpr std::array soil_temp = {"Tsoil[0]"_s, "Tsoil[1]"_s, "Tsoil[2]"_s, "Tsoil[3]"_s};
constexpr int flux_count = 8;
static_assert(sizeof(crocus::Fluxes) == flux_count * sizeof(float));

netCDF::NcVar checkpoint_variable(netCDF::NcFile& file, const std::string& name, std::vector<size_t> shape,
                                  const netCDF::NcType& type)
{
    auto v = file.getVar(name);
    if (v.isNull() || v.getDimCount() != static_cast<int>(shape.size()) || v.getType() != type)
        CHM_THROW_EXCEPTION(module_error, "Missing or incompatible Crocus checkpoint variable " + name);
    for (size_t i = 0; i < shape.size(); ++i)
    {
        if (v.getDim(i).getSize() != shape[i])
            CHM_THROW_EXCEPTION(module_error, "Crocus checkpoint size differs for " + name);
    }
    return v;
}
} // namespace

Crocus::Crocus(config_file cfg)
    : module_base("Crocus", parallel::data, cfg), layers_(cfg.get("max_layers", 20)),
      ground_{cfg.get("ground_temperature_k", 263.f), cfg.get("ground_conductivity", 1.5f),
              cfg.get("ground_layer_thickness_m", .1f), cfg.get("ground_albedo", .2f)}
{
    if (layers_ < 3 ||
        !(ground_.ground_k > 0 && ground_.conductivity > 0 && ground_.ground_layer_m > 0 &&
          ground_.ground_albedo >= 0 && ground_.ground_albedo <= 1) ||
        !std::isfinite(ground_.ground_k) || !std::isfinite(ground_.conductivity) ||
        !std::isfinite(ground_.ground_layer_m))
        CHM_THROW_EXCEPTION(module_error, "Crocus needs at least three layers and a valid ground boundary");

    for (const auto* name :
         {"solar_el", "ilwr", "rh", "t", "p_snow", "p_rain", "U_2m_above_srf", "iswr_direct", "iswr_diffuse"})
        depends(name);
    for (const auto* name :
         {"iswr_subcanopy", "ilwr_subcanopy", "rh_subcanopy", "ta_subcanopy", "p_snow_subcanopy", "p_rain_subcanopy"})
        optional(name);
    for (const auto* name : {"swe",
                             "snowdepthavg",
                             "snowdepthavg_vert",
                             "H",
                             "E",
                             "ilwr_out",
                             "LWout",
                             "melt_rate",
                             "roff",
                             "subl",
                             "sum_snowpack_subl",
                             "snow_albedo",
                             "Nsnow",
                             "Sliq[0]",
                             "Sliq[1]",
                             "Sliq[2]",
                             "Tsnow[0]",
                             "Tsnow[1]",
                             "Tsnow[2]",
                             "Tsoil[0]",
                             "Tsoil[1]",
                             "Tsoil[2]",
                             "Tsoil[3]",
                             "crocus_ground_heat_flux"})
        provides(name);
}

void Crocus::init(mesh& domain)
{
    static std::once_flag constants;
    std::call_once(constants, crocus::chm_crocus_initialize);
    const size_t count = domain->size_local_faces();
    state_.assign(count * crocus::field_count * layers_, 0.f);
    pending_.assign(count, 0.f);
    cumulative_subl_.assign(count, 0.);
    crocus::Fluxes bare{};
    bare.albedo = ground_.ground_albedo;
    bare.lw_out = 5.6703744e-8f * std::pow(ground_.ground_k, 4);
    fluxes_.assign(count, bare);
    pruned_rows_.clear();
    for (size_t i = 0; i < count; ++i)
    {
        auto face = domain->face(i);
        // Point mode retains the original mesh IDs when it prunes unselected faces.
        if (face->cell_local_id != i)
            pruned_rows_.emplace(face->cell_local_id, i);
        publish(face);
    }
}

size_t Crocus::row(const mesh_elem& face) const
{
    const size_t id = face->cell_local_id;
    if (pruned_rows_.empty())
        return id;
    const auto found = pruned_rows_.find(id);
    return found == pruned_rows_.end() ? id : found->second;
}

void Crocus::run(mesh_elem& face)
{
    if (is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }
    const auto get = [&](const char* sub, auto sub_key, auto open_key)
    { return (*face)[has_optional(sub) ? sub_key : open_key]; };
    const double ta = get("ta_subcanopy", "ta_subcanopy"_s, "t"_s);
    const double rh = get("rh_subcanopy", "rh_subcanopy"_s, "rh"_s);
    const double pressure = mio::Atmosphere::stdAirPressure(face->get_z());
    const double vapor = std::clamp(rh, 0., 100.) * .01 * 611.2 * std::exp(17.67 * ta / (ta + 243.5));
    const size_t i = row(face);
    crocus::Forcing m{};
    m.year = global_param->year();
    m.month = global_param->month();
    m.day = global_param->day();
    m.seconds = global_param->hour() * 3600 + global_param->min() * 60 + global_param->sec();
    m.dt = global_param->dt();
    m.air_k = ta + 273.15;
    m.pressure_pa = pressure;
    m.specific_humidity = .622 * vapor / (pressure - .378 * vapor);
    m.wind = std::max(.1, (*face)["U_2m_above_srf"_s]);
    m.sw = has_optional("iswr_subcanopy")
               ? std::max(0., (*face)["iswr_subcanopy"_s])
               : std::max(0., (*face)["iswr_direct"_s]) + std::max(0., (*face)["iswr_diffuse"_s]);
    m.lw = get("ilwr_subcanopy", "ilwr_subcanopy"_s, "ilwr"_s);
    m.snowfall = std::max(0., get("p_snow_subcanopy", "p_snow_subcanopy"_s, "p_snow"_s)) / m.dt;
    m.rainfall = std::max(0., get("p_rain_subcanopy", "p_rain_subcanopy"_s, "p_rain"_s)) / m.dt;
    m.cos_slope = std::max(.001, std::cos(face->slope()));
    m.solar_zenith = (90. - std::clamp((*face)["solar_el"_s], -90., 90.)) * std::numbers::pi / 180.;
    auto* column = state_.data() + i * crocus::field_count * layers_;
    crocus::chm_crocus_step(layers_, column, &m, &ground_, &pending_[i], &fluxes_[i]);
    cumulative_subl_[i] += fluxes_[i].sublimation * m.dt;
    publish(face);
}

void Crocus::publish(mesh_elem& face)
{
    if (is_water(face))
    {
        set_all_nan_on_skip(face);
        return;
    }
    const size_t i = row(face);
    const auto* column = state_.data() + i * crocus::field_count * layers_;
    const auto& flux = fluxes_[i];
    double swe = pending_[i];
    double depth = 0;
    int active = 0;
    while (active < layers_ && column[active] > 0)
    {
        swe += column[active];
        depth += column[active] / column[crocus::density * layers_ + active];
        ++active;
    }
    (*face)["swe"_s] = swe;
    (*face)["snowdepthavg"_s] = depth;
    (*face)["snowdepthavg_vert"_s] = depth / std::max(.001, std::cos(face->slope()));
    (*face)["Nsnow"_s] = active;
    (*face)["H"_s] = flux.sensible;
    (*face)["E"_s] = flux.latent;
    (*face)["LWout"_s] = (*face)["ilwr_out"_s] = flux.lw_out;
    (*face)["melt_rate"_s] = flux.surface_melt;
    (*face)["roff"_s] = flux.runoff;
    (*face)["subl"_s] = flux.sublimation;
    (*face)["sum_snowpack_subl"_s] = cumulative_subl_[i];
    (*face)["snow_albedo"_s] = flux.albedo;
    (*face)["crocus_ground_heat_flux"_s] = flux.basal_heat;
    for (int k = 0; k < 3; ++k)
    {
        // Native H is relative to liquid water at 273.16 K, in J/m2.
        const float specific = k < active ? column[crocus::heat * layers_ + k] / column[k] + 333700.f : 0.f;
        (*face)[snow_temp[k]] = k < active ? 273.16f + std::min(0.f, specific) / 2106.f : -9999.f;
        (*face)[snow_liquid[k]] = k < active ? std::clamp(specific / 333700.f, 0.f, 1.f) * column[k] : 0.f;
    }
    for (auto name : soil_temp)
        (*face)[name] = ground_.ground_k;
}

void Crocus::checkpoint(mesh& domain, netcdf& checkpoint)
{
    auto& file = checkpoint.get_ncfile();
    const auto faces = file.addDim("crocus_face", domain->size_local_faces());
    const auto fields = file.addDim("crocus_field", crocus::field_count);
    const auto layers = file.addDim("crocus_layer", layers_);
    const auto fluxes = file.addDim("crocus_flux", flux_count);
    auto state = file.addVar("crocus_state", netCDF::ncFloat, {faces, fields, layers});
    state.putAtt("schema", netCDF::ncInt, 1);
    state.putAtt("ground_temperature_k", netCDF::ncFloat, ground_.ground_k);
    state.putAtt("fields", "mass density heat diameter sphericity history age");
    state.putVar(state_.data());
    file.addVar("crocus_pending_snow", netCDF::ncFloat, faces).putVar(pending_.data());
    file.addVar("crocus_cumulative_subl", netCDF::ncDouble, faces).putVar(cumulative_subl_.data());
    auto diag = file.addVar("crocus_fluxes", netCDF::ncFloat, {faces, fluxes});
    diag.putAtt("fields", "albedo sensible latent lw_out runoff sublimation surface_melt basal_heat");
    diag.putVar(static_cast<const void*>(fluxes_.data()));
}

void Crocus::load_checkpoint(mesh& domain, netcdf& checkpoint)
{
    auto& file = checkpoint.get_ncfile();
    const auto n = domain->size_local_faces();
    auto state = checkpoint_variable(file, "crocus_state", {n, crocus::field_count, static_cast<size_t>(layers_)},
                                     netCDF::ncFloat);
    int schema;
    float temperature;
    state.getAtt("schema").getValues(&schema);
    state.getAtt("ground_temperature_k").getValues(&temperature);
    if (schema != 1 || temperature != ground_.ground_k)
        CHM_THROW_EXCEPTION(module_error, "Crocus checkpoint schema or fixed ground temperature differs from this run");
    state.getVar(state_.data());
    checkpoint_variable(file, "crocus_pending_snow", {n}, netCDF::ncFloat).getVar(pending_.data());
    checkpoint_variable(file, "crocus_cumulative_subl", {n}, netCDF::ncDouble).getVar(cumulative_subl_.data());
    checkpoint_variable(file, "crocus_fluxes", {n, flux_count}, netCDF::ncFloat)
        .getVar(static_cast<void*>(fluxes_.data()));
    for (size_t i = 0; i < n; ++i)
    {
        auto face = domain->face(i);
        publish(face);
    }
}
