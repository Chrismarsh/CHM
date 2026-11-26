#include "crocus_bridge.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <future>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using namespace crocus;
constexpr int layers = 20;
constexpr Boundary ground{263.f, 1.5f, .1f, .2f};

struct Point
{
    std::array<float, layers * field_count> state{};
    float pending = 0;
    Fluxes flux{ground.ground_albedo};
    float total() const { return std::accumulate(state.begin(), state.begin() + layers, pending); }
    void step(const Forcing& met, const Boundary& b = ground)
    {
        chm_crocus_step(layers, state.data(), &met, &b, &pending, &flux);
    }
};

void near(double got, double expected, double atol, double rtol, const std::string& name)
{
    if (!(std::abs(got - expected) <= atol + rtol * std::abs(expected)))
    {
        std::ostringstream message;
        message << name << " expected " << expected << ", got " << got;
        throw std::runtime_error(message.str());
    }
}

void reference(const std::string& directory)
{
    std::ifstream forcing(directory + "/col_de_porte_forcing.txt");
    std::ifstream reference(directory + "/col_de_porte_reference.txt");
    Forcing met{1989, 12, 20, 0, 3600};
    met.cos_slope = 1;
    met.solar_zenith = 1;
    Point p;
    int step = 0;
    while (forcing >> met.pressure_pa >> met.air_k >> met.specific_humidity >> met.wind >> met.sw >> met.lw >>
           met.snowfall >> met.rainfall)
    {
        met.seconds = ((3 + step) % 24) * 3600;
        met.day = 20 + (3 + step) / 24;
        p.step(met);
        std::string line;
        std::getline(reference, line);
        std::istringstream row(line);
        auto check = [&](double value, double atol, double rtol, const char* name)
        {
            double expected;
            if (!(row >> expected))
                throw std::runtime_error("Missing reference value");
            near(value, expected, atol, rtol, std::to_string(step) + " " + name);
        };
        double depth = 0;
        for (int k = 0; k < layers && p.state[k] > 0; ++k)
            depth += p.state[k] / p.state[density * layers + k];
        check(p.total(), 2e-4, 2e-5, "SWE");
        check(depth, 2e-6, 2e-5, "depth");
        check(p.flux.albedo, 2e-5, 2e-5, "albedo");
        check(p.flux.sensible, 3e-3, 1e-4, "H");
        check(p.flux.runoff, 2e-8, 1e-4, "runoff");
        check(p.flux.sublimation, 2e-9, 1e-4, "sublimation");
        check(p.flux.basal_heat, 3e-3, 1e-4, "ground flux");
        for (int k = 0; k < layers && p.state[k] > 0; ++k)
        {
            const auto value = [&](Field f) { return p.state[f * layers + k]; };
            const float enthalpy = value(heat) / value(mass) + 333700.f;
            check(value(mass), 2e-4, 1e-4, "layer mass");
            check(value(density), 5e-3, 1e-4, "density");
            check(273.16f + std::min(0.f, enthalpy) / 2106.f, .003, 0, "temperature");
            check(std::max(0.f, enthalpy) * value(mass) / 333700000.f, 2e-7, 1e-4, "liquid depth");
            check(value(diameter), 1e-8, 2e-4, "optical diameter");
            check(value(sphericity), 3e-5, 1e-4, "sphericity");
            check(value(history), 1e-5, 0, "grain history");
            check(value(age), 2e-5, 1e-4, "age");
        }
        double extra;
        if (row >> extra)
            throw std::runtime_error("Missing active layer");
        ++step;
    }
    near(step, 24, 0, 0, "number of hours");
}

void water_balance()
{
    Point p;
    Forcing met{2020, 1, 1, 0, 3600, 265, 85000, .002f, 2, 0, 260, 0, 0, 1, 1};
    met.snowfall = .00005f / met.dt;
    met.rainfall = 2.f / met.dt;
    p.step(met);
    near(p.pending, .00005, 1e-10, 0, "trace snowfall");
    near(p.flux.runoff * met.dt, 2, 1e-6, 0, "bare-ground rainfall");
    bool melted = false;
    for (int i = 0; i < 80; ++i)
    {
        const float initial = p.total();
        met.snowfall = i < 10 ? 1.f / met.dt : 0;
        met.rainfall = i == 20 ? 1.f / met.dt : 0;
        if (i >= 20)
        {
            met.air_k = 282;
            met.sw = 500;
            met.lw = 360;
        }
        if (i >= 45)
        {
            met.sw = 950;
            met.lw = 420;
            met.wind = 6;
        }
        p.step(met);
        near(p.total(), initial + (met.snowfall + met.rainfall - p.flux.runoff - p.flux.sublimation) * met.dt, 2e-4,
             2e-5, "column water balance, hour " + std::to_string(i));
        if (!(p.flux.surface_melt >= 0 && std::isfinite(p.flux.lw_out)))
            throw std::runtime_error("Invalid surface flux");
        melted |= p.flux.surface_melt > 0;
    }
    near(p.total(), 0, 0, 0, "complete melt");
    near(melted, true, 0, 0, "melt diagnostic");
    Point warm, cold;
    met.snowfall = 10.f / met.dt;
    met.air_k = 265;
    met.sw = 0;
    met.lw = 260;
    auto warmer = ground;
    warmer.ground_k += 5;
    cold.step(met);
    warm.step(met, warmer);
    if (!(cold.flux.basal_heat > warm.flux.basal_heat))
        throw std::runtime_error("Snow-to-ground heat flux must respond to the prescribed temperature");
}
} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc != 2)
            return 2;
        crocus::chm_crocus_initialize();
        auto other = std::async(std::launch::async, reference, std::string(argv[1]));
        reference(argv[1]);
        water_balance();
        other.get();
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
