#pragma once

namespace crocus
{
// Fortran layout (1, layers, field). Index zero is the top of the snowpack.
enum Field
{
    mass,
    density,
    heat,
    diameter,
    sphericity,
    history,
    age,
    field_count
};
struct Forcing
{
    int year, month, day, seconds;
    float dt, air_k, pressure_pa, specific_humidity, wind, sw, lw;
    float snowfall, rainfall;                    // kg m-2 s-1 at the snow surface
    float cos_slope, solar_zenith;               // radians for zenith
    float transported_mass, transported_density; // positive kg/m2 per step, kg/m3
};

struct Boundary
{
    float ground_k, conductivity, ground_layer_m, ground_albedo;
};

struct Fluxes
{
    // Albedo is prognostic. Supply the previous value on every call.
    float albedo, sensible, latent, lw_out, runoff, sublimation, surface_melt, basal_heat;
};

extern "C" void chm_crocus_initialize();
extern "C" void chm_crocus_step(int layers, float* state, const Forcing* met, const Boundary* ground,
                                float* pending_snow, Fluxes* fluxes);
} // namespace crocus
