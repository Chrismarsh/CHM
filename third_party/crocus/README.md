# ECCC Crocus column

The Fortran sources in `upstream/` come from [ECCC SPS 6.3 at 89a760a](https://github.com/ECCC-ASTD-MRD/sps/tree/89a760a73e860d4972cd61662fb9b6670c657640/src/rpnphy/src/surface/from_surfex). The original copyright headers and CeCILL-C licence are included. CHM builds the externalized (`CROCUS_EXT`, `ECCCGEM`) column directly. The small host module `sfc_options.f90` supplies the RIL exchange option used by the constants initializer.

The C bridge follows the [open-column SNOW_SVS2 call](https://github.com/ECCC-ASTD-MRD/sps/blob/89a760a73e860d4972cd61662fb9b6670c657640/src/rpnphy/src/surface/snow_svs2.F90#L1219). Initialize once before workers start. Each caller owns seven single-precision layer fields, last-step broadband albedo, and a trace snowfall accumulator. Each step supplies a fixed ground temperature, upper ground layer thickness and thermal conductivity. The snow solver computes the resulting heat flux. Wind and temperature reference heights are 2 m. The host supplies terrain-corrected shortwave and precipitation at the snow surface. Crocus resolves metamorphism (B21), albedo/radiation (B92), new snow (V12), thermal conductivity (Y81), retention and compaction (B92). Blowing-snow sublimation is disabled inside this column.

CHM carries the following changes in `snowcro.F90` and its explicit interface in `modi_snowcro.F90`.

- Declare snow albedo `INOUT`. The B92 radiation scheme consumes the previous broadband albedo before the end-of-step update.
- Expose optional surface ice melt in kg m⁻² s⁻¹, measured across the surface melt routine. On complete melt, include remaining ice after sublimation. Drainage separately includes rain, retained meltwater released later, and throughfall on bare ground.
- Route liquid removed above the final holding capacity to drainage. The upstream code computed the excess again after subtracting it, which dropped water during melt. This corrects the flux and preserves the upstream layer update.
- Limit melt to the ice above the single-precision thickness reserve. A layer already below that reserve can otherwise develop negative liquid water during melt and an unrealistic temperature during refreezing.

On a fatal native error, `ABOR1_SFX` flushes the Fortran diagnostics and exits the process with failure so the MPI launcher can terminate the job. This avoids a collective HDF5 file close from an individual rank's exit handlers.

The 24-hour regression starts at hour 03 on model day 354 of 1989 (rows 1899–1922 of [the Col de Porte MESH-SVS forcing](https://github.com/VVionnet/MESH_SVS/blob/aa097436bb2c0981304181b74a1513e98715a035/test_case/basin_forcing.met)). Columns are pressure (Pa), air temperature (K), specific humidity, 2 m wind (m/s, floor 0.1), downwelling shortwave and longwave (W/m²), snowfall and rainfall (kg/m²/s). The reference was generated with **unmodified** SPS `SNOW_SVS2`, 20 layers, an initially bare column, ground reset to 263 K each hour, conductivity 1.5 W/m/K, ground-layer thickness 0.1 m, bare albedo 0.2, flat terrain and zenith 1 rad. Each reference row contains SWE, depth, albedo, sensible heat, drainage, atmospheric mass flux and ground heat, followed by each active layer's SWE, density, temperature, liquid depth, optical diameter, sphericity, history and age. A separate melt and rain test exercises the corrected drainage.

Run the column tests without installing CHM dependencies:

```sh
cmake -S third_party/crocus -B build-crocus
cmake --build build-crocus
ctest --test-dir build-crocus --output-on-failure
```
